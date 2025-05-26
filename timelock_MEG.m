function timelock_MEG(data, save_path, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
timelocked = cell(length(params.trigger_code),1);

% Downsample
if params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    data = ft_resampledata(cfg, data);
end

% Remove padding
cfg = [];
cfg.channel = params.chs;
cfg.latency = [-params.pre params.post];
data = ft_selectdata(cfg, data);

% Demean
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-params.pre 0];
data = ft_preprocessing(cfg,data);

cfg = [];
cfg.demean = 'yes'; %% demean entire trial for whole trial cov
data2 = ft_preprocessing(cfg,data);

for i_phalange = 1:length(params.trigger_code)
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = [-params.pre 0];
    cfg.trials = find(data.trialinfo==params.trigger_code(i_phalange));
    timelocked{i_phalange} = ft_timelockanalysis(cfg, data);
    
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = 'all';
    cfg.trials = find(data2.trialinfo==params.trigger_code(i_phalange));
    tmp = ft_timelockanalysis(cfg, data2);
    timelocked{i_phalange}.cov_all = tmp.cov;
    clear tmp

    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = [0 params.post];
    cfg.trials = find(data.trialinfo==params.trigger_code(i_phalange));
    tmp = ft_timelockanalysis(cfg, data);
    timelocked{i_phalange}.sourcecov = tmp.cov;
    clear tmp

end
save(fullfile(save_path, [params.sub '_' params.modality '_timelocked']), 'timelocked', '-v7.3'); 

%% Find peaks
for i_peak = 1:length(params.peaks)
    peak = cell(length(params.trigger_code),1);
    for i_phalange = 1:length(params.trigger_code)
        dat = timelocked{i_phalange};
        [~, peak_interval(1)] = min(abs(dat.time-params.peaks{i_peak}.peak_latency(1))); % find closest time sample
        [~, peak_interval(2)] = min(abs(dat.time-params.peaks{i_peak}.peak_latency(2))); % find closest time sample
        [~, peak_interval(3)] = min(abs(dat.time-0)); % find closest time sample
        tmp = [];
        t_int = peak_interval(1):peak_interval(2);
        [~, i_peak_latency] = findpeaks(max(dat.avg(:,t_int),[],1)-min(dat.avg(:,t_int),[],1),'SortStr','descend');
        if isempty(i_peak_latency)
            i_peak_latency = round((peak_interval(2)-peak_interval(1))/2);
            tmp.nopeak = true;
        end
        i_peak_latency = peak_interval(1)-1+i_peak_latency(1); % adjust for interval and pick first (=strongest) peak
        tmp.peak_latency = dat.time(i_peak_latency);
        [tmp.max_amplitude, i_maxch] = max(dat.avg(:,i_peak_latency));
        [tmp.min_amplitude, i_minch] = min(dat.avg(:,i_peak_latency));
        tmp.max_channel = dat.label{i_maxch};
        tmp.min_channel = dat.label{i_minch};
        if abs(tmp.max_amplitude) > abs(tmp.min_amplitude)
            tmp.peak_channel = tmp.max_channel;
            tmp.peak_amplitude = abs(tmp.max_amplitude);
            tmp.i_peakch = i_maxch;
        else
            tmp.peak_channel = tmp.min_channel;
            tmp.peak_amplitude = abs(tmp.min_amplitude);
            tmp.i_peakch = i_minch;
        end
        tmp.prestim_std = std(dat.avg(tmp.i_peakch,1:peak_interval(3)));
        tmp.std_error = sqrt(dat.var(tmp.i_peakch,i_peak_latency));
        
        tmp.label = params.peaks{i_peak}.label;
        peak{i_phalange} = tmp;

        %% Plot max channel with variation and peak time
        h = figure;
        hold on
        for i_trl = find(data.trialinfo==params.trigger_code(i_phalange))'
            plot(data.time{i_trl}*1e3, data.trial{i_trl}(peak{i_phalange}.i_peakch,:)*params.amp_scaler,'Color',[211 211 211]/255)
        end
        plot(timelocked{i_phalange}.time*1e3, timelocked{i_phalange}.avg(peak{i_phalange}.i_peakch,:)*params.amp_scaler,'Color',[0 0 0]/255)
        ylimits = ylim;
        lat = 1e3*peak{i_phalange}.peak_latency;
        plot([lat lat],ylimits,'r--')
        hold off
        title(['Peak ' params.modality ' channel: ' peak{i_phalange}.peak_channel])
        ylabel(params.amp_label)
        xlabel('time [ms]')
        xlim([-params.pre params.post]*1e3);
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.peaks{i_peak}.label '_evoked_peakchannel_ph-' params.phalange_labels{i_phalange} '.jpg']))
        close all

        %% Plot topography
        cfg = [];
        cfg.xlim = [peak{i_phalange}.peak_latency-0.005 peak{i_phalange}.peak_latency+0.005];
        cfg.layout = params.layout; 
        cfg.parameter = 'avg';
        h = figure;
        ft_topoplotER(cfg, timelocked{i_phalange});
        axis on
        colorbar
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.peaks{i_peak}.label '_topo_ph-' params.phalange_labels{i_phalange} '.jpg']))
        close all

        %% Plot butterfly
        h = figure;
        plot(timelocked{i_phalange}.time*1e3,timelocked{i_phalange}.avg*params.amp_scaler)
        hold on
        ylimits = ylim;
        lat = 1e3*peak{i_phalange}.peak_latency;
        plot([lat lat],ylimits,'k--')
        hold off
        xlabel('t [msec]')
        ylabel(params.amp_label)
        xlim([-params.pre params.post]*1e3);
        title(['Evoked ' params.modality ' - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(timelocked{i_phalange}.cfg.trials)) ')'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.peaks{i_peak}.label '_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
        close all

        %% Save 
        save(fullfile(save_path, [params.sub '_' params.modality '_' params.peaks{i_peak}.label]), 'peak', '-v7.3'); 
    end
end

end