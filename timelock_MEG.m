function timelock_MEG(data, save_path, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
timelocked = cell(length(params.trigger_codes),1);

% Remove padding
cfg = [];
cfg.toilim = [-params.pre params.post];
data = ft_redefinetrial(cfg, data);

if isfield(params,'baseline')
    baseline = params.baseline;
else
    baseline = [-params.pre 0];
end

% Demean
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = baseline;
data = ft_preprocessing(cfg,data);

if contains(params.chs,'meg')
    cfg = [];
    cfg.channel = 'meg';
    data = ft_selectdata(cfg,data);
else
    cfg = [];
    cfg.channel = params.chs;
    data = ft_selectdata(cfg,data);
end

peak = cell(length(params.trigger_codes),length(params.peaks));
for i_trigger = 1:length(params.trigger_codes)
    if isnumeric(params.trigger_codes{i_trigger}) && size(params.trigger_codes{i_trigger},1)==1 % trigger code(s)
        trls = find(ismember(data.trialinfo,params.trigger_codes{i_trigger}));
        if isfield(params,'maxtrls') && length(trls) > params.maxtrls
            trls = trls(1:params.maxtrls);
        end
    elseif isnumeric(params.trigger_codes{i_trigger}) && size(params.trigger_codes{i_trigger},1)>1 && size(params.trigger_codes{i_trigger},2)==1 % list of trials
        trls = params.trigger_codes{i_trigger};
        if isfield(params,'maxtrls') && length(trls) > params.maxtrls
            trls = trls(1:param.maxtrls);
        end
    elseif ischar(params.trigger_codes{i_trigger}) && strcmp(params.trigger_codes{i_trigger},'all') %
        trls = 1:length(data.trial);
    end
    cfg = [];
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = baseline;
    cfg.trials = trls;
    timelocked{i_trigger} = ft_timelockanalysis(cfg, data);

    vars = mean(timelocked{i_trigger}.var(:,timelocked{i_trigger}.time<=0),2);
    n_sens = length(vars);
    thresh = mean(vars)+3*std(vars);
    h = figure;
    plot(vars);
    hold on
    plot([1 n_sens],[thresh thresh],'k--')
    plot(find(vars>thresh), vars(vars>thresh),'rx')
    text(find(vars>thresh), vars(vars>thresh),timelocked{i_trigger}.label(vars>thresh))
    hold off
    xlabel('Sensor')
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_variances_trig-' params.trigger_labels{i_trigger} '.jpg']))
    close all

    %% Find peaks
    cfg = [];
    cfg.channel = params.chs;
    dat = ft_selectdata(cfg,timelocked{i_trigger});
    if contains(params.chs,'meg')
            cfg = [];
            cfg.channel = params.chs;
            tmp2 = ft_selectdata(cfg,data);
    else
        tmp2 = data;
    end

    for i_peak = 1:length(params.peaks)
        %% Find peak latency and channel
        cfg = [];
        cfg.channel = params.chs;
        dat = ft_selectdata(cfg,timelocked{i_trigger});
        [~, peak_interval(1)] = min(abs(dat.time-params.peaks{i_peak}.peak_latency(1))); % find closest time sample
        [~, peak_interval(2)] = min(abs(dat.time-params.peaks{i_peak}.peak_latency(2))); % find closest time sample
        [~, peak_interval(3)] = min(abs(dat.time-0)); % find closest time sample
        tmp = [];
        t_int = peak_interval(1):peak_interval(2);
        [~, i_peak_latency] = findpeaks(std(dat.avg(:,t_int),0,1),'SortStr','descend');
        if isempty(i_peak_latency)
            lat = (params.peaks{i_peak}.peak_latency(2)-params.peaks{i_peak}.peak_latency(1))/2;
            [~, i_peak_latency] = min(abs(dat.time-lat)); % find closest time sample to center of window
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
        peak{i_trigger,i_peak} = tmp;

        %% Plot max channel with variation and peak time
        h = figure;
        hold on
        for i_trl = trls'
            plot(tmp2.time{i_trl}*1e3, tmp2.trial{i_trl}(peak{i_trigger,i_peak}.i_peakch,:)*params.amp_scaler,'Color',[211 211 211]/255)
        end
        plot(dat.time*1e3, dat.avg(peak{i_trigger,i_peak}.i_peakch,:)*params.amp_scaler,'Color',[0 0 0]/255)
        ylimits = ylim;
        lat = 1e3*peak{i_trigger,i_peak}.peak_latency;
        plot([lat lat],ylimits,'r--')
        hold off
        title(['Peak ' params.modality ' channel: ' peak{i_trigger,i_peak}.peak_channel])
        ylabel(params.amp_label)
        xlabel('time [ms]')
        xlim([-params.pre params.post]*1e3);
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.peaks{i_peak}.label '_evoked_peakchannel_trig-' params.trigger_labels{i_trigger} '.jpg']))
        close all

        %% Plot topography
        cfg = [];
        cfg.xlim = [peak{i_trigger,i_peak}.peak_latency-0.002 peak{i_trigger,i_peak}.peak_latency+0.002];
        cfg.layout = params.layout; 
        cfg.parameter = 'avg';
        cfg.channel = params.chs;
        h = figure;
        ft_topoplotER(cfg, timelocked{i_trigger});
        axis on
        colorbar
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.peaks{i_peak}.label '_topo_trig-' params.trigger_labels{i_trigger} '.jpg']))
        close all
    end
    clear tmp2

    %% Plot butterfly
    h = figure;
    left = 0.1;
    bottom = 0.1;
    width = 0.8;
    gap = 0.10;  % gap between plots
    height_total = 0.8 - gap;  

    height_top = 3/4 * height_total;
    height_bottom = 1/4 * height_total;
    
    % Butterfly
    ax1 = axes('Position', [left, bottom + height_bottom + gap, width, height_top]);
    plot(dat.time*1e3,dat.avg*params.amp_scaler)
    hold on
    ylimits = ylim;
    lat = 1e3*peak{i_trigger,i_peak}.peak_latency;
    plot([lat lat],ylimits,'k--')
    hold off
    xlabel('t [msec]')
    ylabel(params.amp_label)
    xlim([-params.pre params.post]*1e3);
    title(['Evoked ' params.modality ' (n_{trls}=' num2str(length(timelocked{i_trigger}.cfg.trials)) ')'])
    
    % GFP
    ax2 = axes('Position', [left, bottom, width, height_bottom]);
    plot(dat.time*1e3,std(dat.avg,0,1)*params.amp_scaler,'k')
    hold on
    ylimits = ylim;
    lat = 1e3*peak{i_trigger,i_peak}.peak_latency;
    plot([lat lat],ylimits,'k--')
    hold off        
    ax2.XTickLabel = [];
    %xlabel('t [msec]')
    ylabel('GFP')
    xlim([-params.pre params.post]*1e3);
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']))
    close all

    clear dat tmp2
end

%% Save
save(fullfile(save_path, [params.sub '_' params.modality '_timelocked']), 'timelocked', '-v7.3'); 
tmp = peak;
for i_peak = 1:length(params.peaks)
    peak = tmp(:,i_peak);
    save(fullfile(save_path, [params.sub '_' params.modality '_' params.peaks{i_peak}.label]), 'peak', '-v7.3'); 
end

close all

end