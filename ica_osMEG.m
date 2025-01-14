function [data_ica] = ica_osMEG(data,save_path,params)
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
data_ds = ft_resampledata(cfg, data);

%% Calculate ICA
cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = params.n_comp;
cfg.channel    = params.chs;                     
comp           = ft_componentanalysis(cfg, data_ds);

%% Plot
cfg           = [];
cfg.component = 1:n_comp;       
cfg.layout    = params.layout; 
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

savefig(fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_comps'])) 

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
[cfg, ecg_artifact]       = ft_artifact_ecg(cfg,data);

%% Create ECG-locked
cfg = [];
cfg.dftfilter  = 'yes';
cfg.demean     = 'yes';
cfg.trl        = [ecg_artifact zeros(size(ecg_artifact,1), 2)];
temp = ft_redefinetrial(cfg, data);

% Separate MEG and ECG data
cfg.channel    = '*bz';
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
cfg.foilim     = [0 100];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq = ft_freqanalysis(cfg, comp_ecg);
cfg = [];
cfg.channelcmb = {'all' 'ECG'};
cfg.method     = 'coh';
fdcomp = ft_connectivityanalysis(cfg, freq);

% Pick ECG components
maxcoh = max(fdcomp.cohspctrm, [], 2);
ecg_comp_idx = unique(find(R > params.ica_threshold), find(maxcoh > params.ica_threshold));

% Plot correlations
h = figure;
for i = 1:length(ecg_comp_idx)
    subplot(length(ecg_comp_idx),1,i)
    yyaxis left
    plot(timelock.time, timelock.avg(1,:));
    yyaxis right
    plot(timelock.time, timelock.avg(ecg_comp_idx(i)+1,:));  
    title(['Comp: ' num2str(ecg_comp_idx(i)) '; R_{ecg} = ' num2str(R(ecg_comp_idx(i),1))])
end
savefig(fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_ecg_cor'])) 

% Plot coherence spectrum between all components and the ECG
figure;
subplot(3,1,1); plot(fdcomp.freq, abs(fdcomp.cohspctrm)); hold on
plot([min(fdcomp.freq),max(fdcomp.freq)],[param.sica_threshold, parmas.ica_threshold], 'k--')
title('ECG'); xlabel('freq'); ylabel('coh');
subplot(3,1,2); imagesc(abs(fdcomp.cohspctrm));
xlabel('freq'); ylabel('comp');
subplot(3,1,3);
maxcoh = max(fdcomp.cohspctrm, [], 2);
foo = find(~(maxcoh > params.ica_threshold));
bp = bar(1:length(maxcoh), diag(maxcoh), 'stacked');
set(bp(foo),'facecolor','w'); set(bp(ecg_comp_idx),'facecolor','r')
axis([0.5, length(maxcoh)+0.5, 0, 1]); xlabel('comp'); ylabel('coh');

savefig(fullfile(save_path, 'figs',[params.sub '_' params.modality '_ica_ecg_coh'])) 

%% --- EOG ---
% Find EOG artifacts
cfg = [];
cfg.continuous            = 'no';
cfg.channel               = 'EOG';
[~, eog_artifact] = ft_artifact_eog(cfg, data);

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

% Combine EOG and EOG-locked data
comp_eog = ft_appenddata([], eog, comp_eog);

%% Correlation
cfg = [];
timelock = ft_timelockanalysis(cfg, comp_eog);
figure
subplot(2,1,1); plot(timelock.time, timelock.avg(1:2,:)); title('EOG')
subplot(2,1,2); plot(timelock.time, timelock.avg(3:end,:));  title('ICA comp')

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
cfg.foilim     = [0 100];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq = ft_freqanalysis(cfg, comp_eog);
cfg = [];
cfg.method     = 'coh';
cfg.channelcmb = {'comp*' 'EOG001'};
fdcomp_eog1 = ft_connectivityanalysis(cfg, freq);
cfg.channelcmb = {'comp*' 'EOG002'};
fdcomp_eog2 = ft_connectivityanalysis(cfg, freq);

% Find EOG components
maxcoh = max(fdcomp_eog1.cohspctrm, [], 2);
eog1_comp_idx = unique(find(R(:,1) > params.ica_threshold),find(maxcoh > params.ica_threshold));
maxcoh = max(fdcomp_eog2.cohspctrm, [], 2);
eog2_comp_idx = unique(find(R(:,2) > param.sica_threshold),find(maxcoh > params.ica_threshold));

h = figure;
for i = 1:length(eog1_comp_idx)
    subplot(length(eog1_comp_idx),1,i)
    yyaxis left
    plot(timelock.time, timelock.avg(1,:));
    yyaxis right
    plot(timelock.time, timelock.avg(eog1_comp_idx(i)+2,:));  
    title(['Comp: ' num2str(eog1_comp_idx(i)) '; R_{eog1} = ' num2str(R(eog1_comp_idx(i),1))])
end
savefig(fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_eog1_cor'])) 

h = figure;
for i = 1:length(eog2_comp_idx)
    subplot(length(eog2_comp_idx),1,i)
    yyaxis left
    plot(timelock.time, timelock.avg(2,:));
    yyaxis right
    plot(timelock.time, timelock.avg(eog2_comp_idx(i)+2,:));  
    title(['Comp: ' num2str(eog2_comp_idx(i)) '; R_{eog2} = ' num2str(R(eog2_comp_idx(i),2))])
end
savefig(fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_eo2g_cor'])) 

% Plot coherence spectrum between all components and the EOG
figure;
subplot(3,2,1); title('EOG001'); xlabel('freq'); ylabel('coh');
plot(fdcomp_eog1.freq, abs(fdcomp_eog1.cohspctrm)); hold on
plot([min(fdcomp_eog1.freq),max(fdcomp_eog1.freq)],[params.ica_threshold, params.ica_threshold], 'k--');
subplot(3,2,2); title('EOG002'); xlabel('freq'); ylabel('coh');
plot(fdcomp_eog2.freq, abs(fdcomp_eog2.cohspctrm)); hold on
plot([min(fdcomp_eog2.freq),max(fdcomp_eog2.freq)],[params.ica_threshold, params.ica_threshold], 'k--');
subplot(3,2,3); xlabel('freq'); ylabel('comp');
imagesc(abs(fdcomp_eog1.cohspctrm));
subplot(3,2,4); xlabel('freq'); ylabel('comp');
imagesc(abs(fdcomp_eog2.cohspctrm));
subplot(3,2,5); xlabel('comp'); ylabel('coh');
maxcoh = max(fdcomp_eog1.cohspctrm, [], 2);
foo = find(~(maxcoh > params.ica_threshold));
bp = bar(1:length(maxcoh), diag(maxcoh), 'stacked');
set(bp(foo),'facecolor','w'); set(bp(eog1_comp_idx),'facecolor','r');
axis([0.5, length(maxcoh)+0.5, 0, 1]);
subplot(3,2,6); xlabel('comp'); ylabel('coh');
maxcoh = max(fdcomp_eog2.cohspctrm, [], 2);
foo = find(~(maxcoh > params.ica_threshold));
bp = bar(1:length(maxcoh), diag(maxcoh), 'stacked');
set(bp(foo),'facecolor','w'); set(bp(eog2_comp_idx),'facecolor','r'); 
axis([0.5, length(maxcoh)+0.5, 0, 1]);

savefig(fullfile(save_path, 'figs', [params.sub '_' params.modality '_ica_eog_coh'])) 

%% Remove components
% Make a list of all "bad" components
reject_comp = unique([ecg_comp_idx eog1_comp_idx eog2_comp_idx]);

% Remove components
cfg = [];
cfg.component   = reject_comp;
cfg.channel     = params.chs;
cfg.updatesens  = 'no';
data_ica = ft_rejectcomponent(cfg, comp, data);

% Save
save(fullfile(save_path, [params.modality '_ica_comp']), 'comp', 'ecg_comp_idx', 'eog1_comp_idx', 'eog2_comp_idx'); disp('done');
save(fullfile(save_path, [params.modality '_ica']), 'data_ica',"-v7.3"); disp('done');

%%
cfg = [];
cfg.latency = [-0.03 0.3];
data_ica_ds = ft_selectdata(cfg, data_ica);
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-0.03 0];
data_ica_ds = ft_preprocessing(cfg,data_ica_ds);
cfg = [];
cfg.resamplefs = 512;
data_ica_ds = ft_resampledata(cfg, data_ica_ds);
save(fullfile(save_path, [params.modality '_ica_ds']), 'data_ica_ds',"-v7.3"); disp('done');

end