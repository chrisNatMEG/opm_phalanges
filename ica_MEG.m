function [data_ica] = ica_MEG(data, save_path, params, save_results)
%prprocess_osMEG Preprocessing on-scalp MEG data for benchmarking
% recordings. Requires arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), lp_freq, hp_freq,
% bp_freq and notch filter frequencies (corresponding filters are only
% applied if the frequency is defined), n_comp and coh_cutoff (for 
% automated ICA), and ds_freq (downsampling frequency).


%% Downsample
cfg             = [];
cfg.resamplefs  = 200;
cfg.detrend     = 'no';
cfg.channel     = {'EOG', 'ECG', params.chs};
data_ds = ft_resampledata(cfg, data);

%% Calculate ICA
cfg             = [];
cfg.channel  = params.chs;
cfg.trials     = 1;
tmp = ft_selectdata(cfg,data);
data_rank  = rank(tmp.trial{1}*tmp.trial{1}');
clear tmp
if isfield(params,'n_comp') && params.n_comp < data_rank % if n_comp was provided chose the smaller of data_rank and n_comp
   n_comp = params.n_comp;
else
   n_comp = data_rank;
end

cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = n_comp;
cfg.channel    = params.chs;                     
comp           = ft_componentanalysis(cfg, data_ds);

%% Plot all components
if save_results 
    if params.n_comp>16
        for i = 1:floor(n_comp/16)
            cfg           = [];
            cfg.component = (((i-1)*16)+1):(i*16);       
            cfg.layout    = params.layout; 
            cfg.comment   = 'no';
            h = figure;
            ft_topoplotIC(cfg, comp);   
            saveas(h,fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_comps' num2str(i) '.jpg'])) 
        end
    
        if mod(params.n_comp,16)~=0
            cfg           = [];
            cfg.component = ((i*16)+1):n_comp;       
            cfg.layout    = params.layout; 
            cfg.comment   = 'no';
            h = figure;
            ft_topoplotIC(cfg, comp);   
            saveas(h,fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_comps' num2str(i+1) '.jpg'])) 
        end
    else
        cfg           = [];
        cfg.component = 1:n_comp;       
        cfg.layout    = params.layout; 
        cfg.comment   = 'no';
        h = figure;
        ft_topoplotIC(cfg, comp);   
        saveas(h,fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_comps.jpg'])) 
    end
    close all
end

%% --- ECG ----
% Find ECG artifacts
cfg                       = [];
cfg.continuous            = 'no';
cfg.artfctdef.ecg.pretim  = 0.25;
cfg.artfctdef.ecg.psttim  = 0.40-1/1200;
cfg.artfctdef.ecg.cutoff  = 3;
cfg.artfctdef.ecg.feedback = 'no';
cfg.channel               = 'ECG';
cfg.artfctdef.ecg.inspect = 'ECG';
[cfg, ecg_artifact]       = ft_artifact_ecg(cfg,data); % if ft_fetch_data creates errors add "'allowoverlap', true" to the function call.

%% Create ECG-locked
cfg = [];
cfg.dftfilter  = 'yes';
cfg.demean     = 'yes';
cfg.trl        = [ecg_artifact zeros(size(ecg_artifact,1), 2)];
temp = ft_redefinetrial(cfg, data);

% Separate MEG and ECG data
cfg.channel    = params.chs;
data_ecg = ft_selectdata(cfg, temp);
cfg.channel    = 'ECG';
ecg = ft_selectdata(cfg, temp);
ecg.channel{:} = 'ECG';

% Filter power line noise
cfg = [];
cfg.dftfilter       = 'yes';
cfg.dftfreq         = [50, 100, 150];
ecg = ft_preprocessing(cfg, ecg);

% Decompose the ECG-locked data
cfg = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_ecg = ft_componentanalysis(cfg, data_ecg);

% Combine ECG and ECG-locked data
comp_ecg = ft_appenddata([], ecg, comp_ecg);

cfg = [];
timelock = ft_timelockanalysis(cfg, comp_ecg);

%% Correlation
ecg_comp_idx = [];
for i = 2:size(timelock.avg,1)
    tmp = corrcoef(timelock.avg(1,:), timelock.avg(i,:));
    R(i-1,1) = tmp(2,1);
end

% Find components with high ECG coherence
% Compute coherence between all components and the ECG
cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.foilim     = [0 params.filter.lp_freq];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq = ft_freqanalysis(cfg, comp_ecg);
cfg = [];
cfg.channelcmb = {'all' 'ECG'};
cfg.method     = 'coh';
fdcomp = ft_connectivityanalysis(cfg, freq);
clear freq

% Pick ECG components
maxcoh = max(fdcomp.cohspctrm, [], 2);
ecg_comp_idx = unique([find(abs(R) > params.ica_cor); find(maxcoh > params.ica_coh)]);

if save_results
    % Plot correlations
    if length(ecg_comp_idx)>=1
        h = figure;
        for i = 1:length(ecg_comp_idx)
            subplot(length(ecg_comp_idx),1,i)
            yyaxis left
            plot(timelock.time, timelock.avg(1,:));
            yyaxis right
            plot(timelock.time, timelock.avg(ecg_comp_idx(i)+1,:));  
            title(['Comp: ' num2str(ecg_comp_idx(i)) '; R_{ecg} = ' num2str(R(ecg_comp_idx(i),1))])
        end
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_ecg_cor.jpg'])) 
    end
    close all
    
    % Plot coherence spectrum between all components and the ECG
    h = figure;
    subplot(3,1,1); plot(fdcomp.freq, abs(fdcomp.cohspctrm)); hold on
    plot([min(fdcomp.freq),max(fdcomp.freq)],[params.ica_coh, params.ica_coh], 'k--')
    title('ECG'); xlabel('freq'); ylabel('coh');
    subplot(3,1,2); imagesc(abs(fdcomp.cohspctrm));
    xlabel('freq'); ylabel('comp');
    subplot(3,1,3);
    maxcoh = max(fdcomp.cohspctrm, [], 2);
    foo = find(~(maxcoh > params.ica_coh));
    bp = bar(1:length(maxcoh), diag(maxcoh), 'stacked');
    set(bp(foo),'facecolor','w'); set(bp(ecg_comp_idx),'facecolor','r')
    axis([0.5, length(maxcoh)+0.5, 0, 1]); xlabel('comp'); ylabel('coh');    
    saveas(h, fullfile(save_path, 'figs',[params.sub '_' params.modality '_ica_ecg_coh.jpg'])) 
    close all
end

%% --- EOG ---
% Find EOG artifacts
cfg = [];
cfg.continuous            = 'no';
cfg.channel               = 'EOG';
[~, eog_artifact] = ft_artifact_eog(cfg, data); % if ft_fetch_data creates errors add "'allowoverlap', true" to the function call.

% Make artifact epochs
cfg = [];
cfg.dftfilter  = 'yes';
cfg.demean     = 'yes';
cfg.trl        = [eog_artifact zeros(size(eog_artifact,1), 1)];
temp = ft_redefinetrial(cfg, data);
    
% Separate MEG and EOG data
cfg.channel    = params.chs;
data_eog = ft_selectdata(cfg, temp);
cfg.channel    = 'EOG';
eog = ft_selectdata(cfg, temp);
eog.channel{:} = 'EOG';         % renaming for bookkeeping

% Exclude too short trials
for i = 1:length(eog.trial)
    incl(i) = size(eog.trial{i},2)>10;
end
cfg = [];
cfg.trials = incl;
eog = ft_selectdata(cfg,eog);

% Filter power line noise
cfg = [];
cfg.dftfilter  = 'yes';
cfg.dftfreq    = [50, 100, 150];
eog = ft_preprocessing(cfg, eog);

% Decompose EOG-locked data
cfg = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_eog = ft_componentanalysis(cfg, data_eog);

cfg = [];
cfg.trials = incl;
comp_eog = ft_selectdata(cfg,comp_eog);

% Combine EOG and EOG-locked data
comp_eog = ft_appenddata([], eog, comp_eog);

%% Correlation
cfg = [];
timelock = ft_timelockanalysis(cfg, comp_eog);

eog1_comp_idx = [];
eog2_comp_idx = [];
for i = 3:size(timelock.avg,1)
    tmp = corrcoef(timelock.avg(1,:), timelock.avg(i,:));
    R(i-2,1) = tmp(2,1);
    tmp = corrcoef(timelock.avg(2,:), timelock.avg(i,:));
    R(i-2,2) = tmp(2,1);
end

%% Compute coherence between all components and the E0G
cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.foilim     = [0 params.filter.lp_freq];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq = ft_freqanalysis(cfg, comp_eog);
cfg = [];
cfg.method     = 'coh';
cfg.channelcmb = {'comp*' 'EOG001'};
fdcomp_eog1 = ft_connectivityanalysis(cfg, freq);
cfg.channelcmb = {'comp*' 'EOG002'};
fdcomp_eog2 = ft_connectivityanalysis(cfg, freq);
clear freq

% Find EOG components
maxcoh = max(fdcomp_eog1.cohspctrm, [], 2);
eog1_comp_idx = unique([find(abs(R(:,1)) > params.ica_cor); find(maxcoh > params.ica_coh)]);
maxcoh = max(fdcomp_eog2.cohspctrm, [], 2);
eog2_comp_idx = unique([find(abs(R(:,2)) > params.ica_cor); find(maxcoh > params.ica_coh)]);

if save_results
    if length(eog1_comp_idx)>=1
        h = figure;
        for i = 1:length(eog1_comp_idx)
            subplot(length(eog1_comp_idx),1,i)
            yyaxis left
            plot(timelock.time, timelock.avg(1,:));
            yyaxis right
            plot(timelock.time, timelock.avg(eog1_comp_idx(i)+2,:));  
            title(['Comp: ' num2str(eog1_comp_idx(i)) '; R_{eog1} = ' num2str(R(eog1_comp_idx(i),1))])
        end
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_eog1_cor.jpg'])) 
    end
    close all

    if length(eog2_comp_idx)>=1
        h = figure;
        for i = 1:length(eog2_comp_idx)
            subplot(length(eog2_comp_idx),1,i)
            yyaxis left
            plot(timelock.time, timelock.avg(2,:));
            yyaxis right
            plot(timelock.time, timelock.avg(eog2_comp_idx(i)+2,:));  
            title(['Comp: ' num2str(eog2_comp_idx(i)) '; R_{eog2} = ' num2str(R(eog2_comp_idx(i),2))])
        end
        saveas(h,fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_eog2_cor.jpg'])) 
    end
    close all

    % Plot coherence spectrum between all components and the EOG
    h = figure;
    subplot(3,2,1); title('EOG001'); xlabel('freq'); ylabel('coh');
    plot(fdcomp_eog1.freq, abs(fdcomp_eog1.cohspctrm)); hold on
    plot([min(fdcomp_eog1.freq),max(fdcomp_eog1.freq)],[params.ica_coh, params.ica_coh], 'k--');
    subplot(3,2,2); title('EOG002'); xlabel('freq'); ylabel('coh');
    plot(fdcomp_eog2.freq, abs(fdcomp_eog2.cohspctrm)); hold on
    plot([min(fdcomp_eog2.freq),max(fdcomp_eog2.freq)],[params.ica_coh, params.ica_coh], 'k--');
    subplot(3,2,3); xlabel('freq'); ylabel('comp');
    imagesc(abs(fdcomp_eog1.cohspctrm));
    subplot(3,2,4); xlabel('freq'); ylabel('comp');
    imagesc(abs(fdcomp_eog2.cohspctrm));
    subplot(3,2,5); xlabel('comp'); ylabel('coh');
    maxcoh = max(fdcomp_eog1.cohspctrm, [], 2);
    foo = find(~(maxcoh > params.ica_coh));
    bp = bar(1:length(maxcoh), diag(maxcoh), 'stacked');
    set(bp(foo),'facecolor','w'); set(bp(eog1_comp_idx),'facecolor','r');
    axis([0.5, length(maxcoh)+0.5, 0, 1]);
    subplot(3,2,6); xlabel('comp'); ylabel('coh');
    maxcoh = max(fdcomp_eog2.cohspctrm, [], 2);
    foo = find(~(maxcoh > params.ica_coh));
    bp = bar(1:length(maxcoh), diag(maxcoh), 'stacked');
    set(bp(foo),'facecolor','w'); set(bp(eog2_comp_idx),'facecolor','r'); 
    axis([0.5, length(maxcoh)+0.5, 0, 1]);
    
    saveas(h,fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_eog_coh.jpg'])) 
    close all
end

%% Remove components
% Make a list of all "bad" components
reject_comp = unique([ecg_comp_idx; eog1_comp_idx; eog2_comp_idx]);

% Remove components
cfg = [];
cfg.component   = reject_comp;
cfg.channel     = params.chs;
cfg.updatesens  = 'no';
data_ica = ft_rejectcomponent(cfg, comp, data);

%% Save components
if save_results
    save(fullfile(save_path, [params.sub '_' params.modality '_ica_comp']), 'comp', 'ecg_comp_idx', 'eog1_comp_idx', 'eog2_comp_idx'); disp('done');
    
    cfg           = [];
    cfg.component = reject_comp;       
    cfg.layout    = params.layout; 
    cfg.comment   = 'no';
    
    if length(reject_comp)>=1
        h = figure;
        ft_topoplotIC(cfg, comp);   
        saveas(h,fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_rejected_comps' num2str(i) '.jpg'])) 
    end

    %% Downsample for MVPA
    cfg = [];
    cfg.latency = [-0.05 0.3];
    data_ica_ds = ft_selectdata(cfg, data_ica);
    cfg = [];
    cfg.demean = 'yes';
    cfg.baselinewindow = [-0.05 0];
    data_ica_ds = ft_preprocessing(cfg,data_ica_ds);
    cfg = [];
    cfg.resamplefs = 200;
    data_ica_ds = ft_resampledata(cfg, data_ica_ds);
    save(fullfile(save_path, [params.sub '_' params.modality '_ica_ds']), 'data_ica_ds',"-v7.3"); disp('done');
    clear data_ica_ds
end

close all
end