function sensor_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n_triggers = length(params.trigger_labels);

for i_peak = 1:length(params.peaks)
    peak_label = ['_' params.peaks{i_peak}.label];
    
    peak_ratio = [];
    snr = [];
    latency = [];
    n_trl = [];
    for i_sub = subs
        params.sub = ['sub_' num2str(i_sub,'%02d')];
        ft_hastoolbox('mne', 1);
        save_path = fullfile(base_save_path,params.sub);
        clear peak_opm peak_opmeeg
        clear peak_squid peak_squideeg
        peak_opm = load(fullfile(save_path, [params.sub '_opm_' params.peaks{1}.label '.mat'])).peak; 
        peak_opmeeg = load(fullfile(save_path, [params.sub '_opmeeg_' params.peaks{1}.label '.mat'])).peak; 
        peak_squid = load(fullfile(save_path, [params.sub '_squid_' params.peaks{1}.label '.mat'])).peak; 
        peak_squideeg = load(fullfile(save_path, [params.sub '_squideeg_' params.peaks{1}.label '.mat'])).peak;
    
        clear squid_timelocked opm_timelocked squideeg_timelocked opmeeg_timelocked
        squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked'])).timelocked; 
        opm_timelocked = load(fullfile(save_path, [params.sub '_opm_timelocked'])).timelocked; 
        squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked'])).timelocked; 
        opmeeg_timelocked = load(fullfile(save_path, [params.sub '_opmeeg_timelocked'])).timelocked; 
        meg_chs = find(contains(squid_timelocked{1}.label,'MEG')&endsWith(squid_timelocked{1}.label,'1'));
        opm_chs = find(contains(opm_timelocked{1}.label,'bz'));
    
        for i_trigger = 1:n_triggers
            peak_ratio.meg(i_sub,i_trigger) = peak_opm{i_trigger}.peak_amplitude/peak_squid{i_trigger}.peak_amplitude;
            peak_ratio.eeg(i_sub,i_trigger) = peak_opmeeg{i_trigger}.peak_amplitude/peak_squideeg{i_trigger}.peak_amplitude;
            snr.error_opm(i_sub,i_trigger) = peak_opm{i_trigger}.peak_amplitude/peak_opm{i_trigger}.std_error;
            snr.error_squidmag(i_sub,i_trigger) = peak_squid{i_trigger}.peak_amplitude/peak_squid{i_trigger}.std_error;
            snr.error_squideeg(i_sub,i_trigger) = peak_squideeg{i_trigger}.peak_amplitude/peak_squideeg{i_trigger}.std_error;
            snr.error_opmeeg(i_sub,i_trigger) = peak_opmeeg{i_trigger}.peak_amplitude/peak_opmeeg{i_trigger}.std_error;
            snr.prestim_opm(i_sub,i_trigger) = peak_opm{i_trigger}.peak_amplitude/peak_opm{i_trigger}.prestim_std;
            snr.prestim_squidmag(i_sub,i_trigger) = peak_squid{i_trigger}.peak_amplitude/peak_squid{i_trigger}.prestim_std;
            snr.prestim_opmeeg(i_sub,i_trigger) = peak_opmeeg{i_trigger}.peak_amplitude/peak_opmeeg{i_trigger}.prestim_std;
            snr.prestim_squideeg(i_sub,i_trigger) = peak_squideeg{i_trigger}.peak_amplitude/peak_squideeg{i_trigger}.prestim_std;
            snr.ratio_error(i_sub,i_trigger) = snr.error_opm(i_sub,i_trigger)/snr.error_squidmag(i_sub,i_trigger);
            snr.ratio_prestim(i_sub,i_trigger) = snr.prestim_opm(i_sub,i_trigger)/snr.prestim_squidmag(i_sub,i_trigger);
            latency.opm(i_sub,i_trigger) = peak_opm{i_trigger}.peak_latency;
            latency.squidmag(i_sub,i_trigger) = peak_squid{i_trigger}.peak_latency;
            latency.opmeeg(i_sub,i_trigger) = peak_opmeeg{i_trigger}.peak_latency;
            latency.squideeg(i_sub,i_trigger) = peak_squideeg{i_trigger}.peak_latency;
            amp.opm(i_sub,i_trigger) = peak_opm{i_trigger}.peak_amplitude;
            amp.squidmag(i_sub,i_trigger) = peak_squid{i_trigger}.peak_amplitude;
            amp.opmeeg(i_sub,i_trigger) = peak_opmeeg{i_trigger}.peak_amplitude;
            amp.squideeg(i_sub,i_trigger) = peak_squideeg{i_trigger}.peak_amplitude;
            n_trl.opm(i_sub,i_trigger) = length(opm_timelocked{i_trigger}.cfg.trials);
            n_trl.squidmag(i_sub,i_trigger) = length(squid_timelocked{i_trigger}.cfg.trials);
            n_trl.opmeeg(i_sub,i_trigger) = length(opmeeg_timelocked{i_trigger}.cfg.trials);
            n_trl.squideeg(i_sub,i_trigger) = length(squideeg_timelocked{i_trigger}.cfg.trials);
    
            h = figure;
            subplot(2,1,1)
            plot(squid_timelocked{i_trigger}.time*1e3,squid_timelocked{i_trigger}.avg(meg_chs(1:3:end),:)*1e15)
            xlabel('t [msec]')
            ylabel('B [fT]')
            xlim([-params.pre params.post]*1e3)
            title(['Evoked SQUID MAG - ' params.trigger_labels{i_trigger} ' (n_{trls}=' num2str(length(squid_timelocked{i_trigger}.cfg.trials)) ')'])
            subplot(2,1,2)
            plot(opm_timelocked{i_trigger}.time*1e3,opm_timelocked{i_trigger}.avg(opm_chs,:)*1e15)
            xlabel('t [msec]')
            ylabel('B [fT]')
            xlim([-params.pre params.post]*1e3)
            title(['Evoked OPM MAG - ' params.trigger_labels{i_trigger} ' (n_{trls}=' num2str(length(opm_timelocked{i_trigger}.cfg.trials)) ')'])
            saveas(h, fullfile(save_path, 'figs', [params.sub '_squidopm_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']))
            close all
        end
    end
    clear squid_timelocked opm_timelocked squideeg_timelocked opmeeg_timelocked
    
    %% Save
    save(fullfile(base_save_path, ['group_sensor' peak_label]),"peak_ratio","snr","latency","amp","n_trl","-v7.3");
    
    %% Plot ratio
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,mean(peak_ratio.meg,1,'omitnan'));
    hold on
    er = errorbar(1:n_triggers,mean(peak_ratio.meg,1,'omitnan'), mean(peak_ratio.meg,1,'omitnan')-min(peak_ratio.meg,[],1,'omitnan'), mean(peak_ratio.meg,1,'omitnan')-max(peak_ratio.meg,[],1,'omitnan'));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    er.LineWidth = 1;
    er.CapSize = 30;
    hold off
    title([params.peaks{1}.label ' peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.meg,'omitnan'),'omitnan'),'%.2f') ')'])
    ylabel('OPM/SQUID')
    xlabel('Trigger')
    xticklabels(params.trigger_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['Peak_amplitude_ratios_meg' peak_label '.jpg']))
    close all
    
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,mean(peak_ratio.eeg,1,'omitnan'));
    hold on
    er = errorbar(1:n_triggers,mean(peak_ratio.eeg,1,'omitnan'), mean(peak_ratio.eeg,1,'omitnan')-min(peak_ratio.eeg,[],1,'omitnan'), mean(peak_ratio.eeg,1,'omitnan')-max(peak_ratio.eeg,[],1,'omitnan'));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    er.LineWidth = 1;
    er.CapSize = 30;
    hold off
    title([params.peaks{1}.label ' peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.eeg,'omitnan'),'omitnan'),'%.2f') ')'])
    ylabel('OPMEEG/SQUIDEEG')
    xlabel('Trigger')
    xticklabels(params.trigger_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['Peak_amplitude_ratios_eeg' peak_label '.jpg']))
    close all
    
    %% Plot SNR - error
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,mean(snr.ratio_error,1,'omitnan'));
    hold on
    er = errorbar(1:n_triggers,mean(snr.ratio_error,1,'omitnan'), mean(snr.ratio_error,1,'omitnan')-min(snr.ratio_error,[],1,'omitnan'), mean(snr.ratio_error,1,'omitnan')-max(snr.ratio_error,[],1,'omitnan'));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    er.LineWidth = 1;
    er.CapSize = 30;
    hold off
    title([params.peaks{1}.label ' SNR_{stderror} ratio (mean = ' num2str(mean(mean(snr.ratio_error,'omitnan'),'omitnan'),'%.2f') ')'])
    ylabel('OPM/SQUID')
    xlabel('Trigger')
    xticklabels(params.trigger_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['SNR_ratios_error' peak_label '.jpg']))
    close all
    
    %% Plot SNR - prestim
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,mean(snr.ratio_prestim,1,'omitnan'));
    hold on
    er = errorbar(1:n_triggers,mean(snr.ratio_prestim,1,'omitnan'), mean(snr.ratio_prestim,1,'omitnan')-min(snr.ratio_prestim,[],1,'omitnan'), mean(snr.ratio_prestim,1,'omitnan')-max(snr.ratio_prestim,[],1,'omitnan'));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    er.LineWidth = 1;
    er.CapSize = 30;
    hold off
    title([params.peaks{1}.label ' SNR_{prestim} ratio (mean = ' num2str(mean(mean(snr.ratio_prestim,'omitnan'),'omitnan'),'%.2f') ')'])
    ylabel('OPM/SQUID')
    xlabel('Trigger')
    xticklabels(params.trigger_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['SNR_ratios_prestim' peak_label '.jpg']))
    close all
    
    %% Plot peak amp
    % MEG
    data1 = 1e15*amp.squidmag;
    data2 = 1e15*amp.opm;
    mean1 = mean(data1,1,'omitnan');
    mean2 = mean(data2,1,'omitnan');
    min1 = min(data1,[],1,'omitnan');
    min2 = min(data2,[],1,'omitnan');
    max1 = max(data1,[],1,'omitnan');
    max2 = max(data2,[],1,'omitnan');
    err1 = [mean1-min1; max1-mean1];
    err2 = [mean2-min2; max2-mean2];
    
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,[mean1; mean2]','grouped');
    hold on
    for k=1:n_triggers
        errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    end

    p_values = zeros(1,n_triggers);
    for i = 1:n_triggers
        [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
    end
    sigstar(arrayfun(@(x) [x, x], 1:n_triggers, 'UniformOutput', false), p_values);
    hold off
    title(['Group level ' params.peaks{1}.label ' amplitude'])
    ylabel('Peak amplitude [fT]')
    xlabel('Trigger')
    xticklabels(params.trigger_labels)
    legend({'squidmag','opm'});
    saveas(h, fullfile(base_save_path, 'figs', ['Amplitude_meg' peak_label '.jpg']))
    close all 
    
    data = {data1, data2};
    triggerLabels = params.trigger_labels;
    yLabelStr = 'Peak amplitude [fT]';
    titleStr = ['Group level ' params.peaks{1}.label ' amplitude'];
    save_path = fullfile(base_save_path, 'figs', ['Amplitude_meg_box' peak_label '.jpg']);
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
    %% EEG
    data1 = 1e6*amp.squideeg;
    data2 = 1e6*amp.opmeeg;
    mean1 = mean(data1,1,'omitnan');
    mean2 = mean(data2,1,'omitnan');
    min1 = min(data1,[],1,'omitnan');
    min2 = min(data2,[],1,'omitnan');
    max1 = max(data1,[],1,'omitnan');
    max2 = max(data2,[],1,'omitnan');
    err1 = [mean1-min1; max1-mean1];
    err2 = [mean2-min2; max2-mean2];
    
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,[mean1; mean2]','grouped');
    hold on
    for k=1:n_triggers
        errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    end
    p_values = zeros(1,n_triggers);
    for i = 1:n_triggers
        [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
    end
    sigstar(arrayfun(@(x) [x, x], 1:n_triggers, 'UniformOutput', false), p_values);
    hold off
    title(['Group level ' params.peaks{1}.label ' amplitude'])
    ylabel('Peak amplitude [uV]')
    xlabel('Trigger')
    xticklabels(params.trigger_labels)
    legend({'squideeg','opmeeg'});
    saveas(h, fullfile(base_save_path, 'figs', ['Amplitude_eeg' peak_label '.jpg']))
    close all
    
    data = {data1, data2};
    triggerLabels = params.trigger_labels;
    yLabelStr = 'Peak amplitude [uV]';
    titleStr = ['Group level ' params.peaks{1}.label ' amplitude'];
    save_path = fullfile(base_save_path, 'figs', ['Amplitude_eeg_box' peak_label '.jpg']);
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
    
    %% Plot peak latency
    data1 = 1e3*latency.squidmag;
    data2 = 1e3*latency.opm;
    mean1 = mean(data1,1,'omitnan');
    mean2 = mean(data2,1,'omitnan');
    min1 = min(data1,[],1,'omitnan');
    min2 = min(data2,[],1,'omitnan');
    max1 = max(data1,[],1,'omitnan');
    max2 = max(data2,[],1,'omitnan');
    err1 = [mean1-min1; max1-mean1];
    err2 = [mean2-min2; max2-mean2];
    
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,[mean1; mean2]','grouped');
    hold on
    for k=1:n_triggers
        errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    end
    p_values = zeros(1,n_triggers);
    for i = 1:n_triggers
        [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
    end
    sigstar(arrayfun(@(x) [x, x], 1:n_triggers, 'UniformOutput', false), p_values);
    hold off
    title(['Group level ' params.peaks{1}.label ' latency'])
    ylabel('Latency [ms]')
    xlabel('Trigger')
    xticklabels(params.trigger_labels)
    legend({'squidmag','opm'},'Location','southeast');
    saveas(h, fullfile(base_save_path, 'figs', ['Latency' peak_label '.jpg']))
    close all
    
    data = {data1, data2};
    triggerLabels = params.trigger_labels;
    yLabelStr = 'Latency [ms]';
    titleStr = ['Group level ' params.peaks{1}.label ' latency'];
    save_path = fullfile(base_save_path, 'figs', ['Latency_box' peak_label '.jpg']);
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
    
    %% Plot SNR - error
    data1 = snr.error_squidmag;
    data2 = snr.error_opm;
    mean1 = mean(data1,1,'omitnan');
    mean2 = mean(data2,1,'omitnan');
    min1 = min(data1,[],1,'omitnan');
    min2 = min(data2,[],1,'omitnan');
    max1 = max(data1,[],1,'omitnan');
    max2 = max(data2,[],1,'omitnan');
    err1 = [mean1-min1; max1-mean1];
    err2 = [mean2-min2; max2-mean2];
    
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,[mean1; mean2]','grouped');
    hold on
    for k=1:n_triggers
        errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    end
    p_values = zeros(1,n_triggers);
    for i = 1:n_triggers
        [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
    end
    sigstar(arrayfun(@(x) [x, x], 1:n_triggers, 'UniformOutput', false), p_values);
    hold off
    title('Group level SNR_{stderror}')
    ylabel('SNR')
    xlabel('Trigger')
    legend({'squidmag','opm'});
    xticklabels(params.trigger_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['SNR_error' peak_label '.jpg']))
    close all
    
    data = {data1, data2};
    triggerLabels = params.trigger_labels;
    yLabelStr = 'SNR';
    titleStr = ['Group level ' params.peaks{1}.label ' SNR_{stderror}'];
    save_path = fullfile(base_save_path, 'figs', ['SNR_error_box' peak_label '.jpg']);
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
    
    %% Plot SNR - prestim
    data1 = snr.prestim_squidmag;
    data2 = snr.prestim_opm;
    mean1 = mean(data1,1,'omitnan');
    mean2 = mean(data2,1,'omitnan');
    min1 = min(data1,[],1,'omitnan');
    min2 = min(data2,[],1,'omitnan');
    max1 = max(data1,[],1,'omitnan');
    max2 = max(data2,[],1,'omitnan');
    err1 = [mean1-min1; max1-mean1];
    err2 = [mean2-min2; max2-mean2];
    
    h = figure('DefaultAxesFontSize',16);
    bar(1:n_triggers,[mean1; mean2]','grouped');
    hold on
    for k=1:n_triggers
        errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    end
    p_values = zeros(1,n_triggers);
    for i = 1:n_triggers
        [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
    end
    sigstar(arrayfun(@(x) [x, x], 1:n_triggers, 'UniformOutput', false), p_values);
    hold off
    title('Group level SNR_{prestim}')
    ylabel('SNR')
    xlabel('Trigger')
    legend({'squidmag','opm'});
    xticklabels(params.trigger_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['SNR_prestim' peak_label '.jpg']))
    close all
    
    data = {data1, data2};
    triggerLabels = params.trigger_labels;
    yLabelStr = 'SNR';
    titleStr = ['Group level ' params.peaks{1}.label ' SNR_{prestim}'];
    save_path = fullfile(base_save_path, 'figs', ['SNR_prestim_box' peak_label '.jpg']);
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path, 1);

    %% Number of trials
    data1 = n_trl.squideeg;
    data2 = n_trl.opmeeg;
    data = {data1, data2};
    triggerLabels = params.trigger_labels;
    yLabelStr = 'n_{trl}';
    titleStr = ['Group level number of trials'];
    save_path = fullfile(base_save_path, 'figs', ['Ntrl_EEG_box' peak_label '.jpg']);
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
    data1 = n_trl.squidmag;
    data2 = n_trl.opm;
    data = {data1, data2};
    triggerLabels = params.trigger_labels;
    yLabelStr = 'n_{trl}';
    titleStr = ['Group level number of trials'];
    save_path = fullfile(base_save_path, 'figs', ['Ntrl_MEG_box' peak_label '.jpg']);
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

end

%% Grand average
squid = cell(n_triggers,length(subs));
opm = cell(n_triggers,length(subs));
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);

    squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked'])).timelocked; 
    opm_timelocked = load(fullfile(save_path, [params.sub '_opm_timelocked'])).timelocked; 

    for i_trigger = 1:n_triggers
        squid{i_trigger,i_sub} = squid_timelocked{i_trigger};
        opm{i_trigger,i_sub} = opm_timelocked{i_trigger};
    end
    clear squid_timelocked opm_timelocked
end
for i_trigger = 1:n_triggers

    % OPM
    cfg = [];
    cfg.channel = '*bz';
    grandavg_opm = ft_timelockgrandaverage(cfg,opm{i_trigger,~cellfun(@isempty,opm(i_trigger,:))});

    h = figure; 
    plot(grandavg_opm.time*1e3,grandavg_opm.avg*1e15)
    xlabel('t [msec]')
    ylabel('B [fT]')
    xlim([-params.pre params.post]*1e3);
    title(['Grand average opm - phalange ' params.trigger_labels{i_trigger}])
    saveas(h, fullfile(base_save_path, 'figs', ['opm_grndAvg_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']))
    close all
    
    for i_peak = 1:length(params.peaks)
        peak_label = ['_' params.peaks{i_peak}.label];
        cfg = [];
        cfg.xlim = params.peaks{1}.peak_latency;
        %cfg.zlim = [0 6e-14];
        cfg.layout = 'fieldlinebeta2bz_helmet.mat';
        h = figure; ft_topoplotER(cfg,grandavg_opm); colorbar; title(['GRAND AVG OPM - ' params.trigger_labels{i_trigger}])
        saveas(h, fullfile(base_save_path, 'figs', ['opm' peak_label '_grndAvg_topo_trig-' params.trigger_labels{i_trigger} '.jpg']))
        close all
    end

    % SQUID-GRAD
    cfg = [];
    cfg.channel = 'meggrad';
    grandavg_squidgrad = ft_timelockgrandaverage(cfg,squid{i_trigger,~cellfun(@isempty,squid(i_trigger,:))});

    h = figure; 
    plot(grandavg_squidgrad.time*1e3,grandavg_squidgrad.avg*1e15)
    xlabel('t [msec]')
    ylabel('B [fT]')
    xlim([-params.pre params.post]*1e3);
    title(['Grand average squidgrad - phalange ' params.trigger_labels{i_trigger}])
    saveas(h, fullfile(base_save_path, 'figs', ['squidgrad_grndAvg_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']))
    close all

    for i_peak = 1:length(params.peaks)
        peak_label = ['_' params.peaks{i_peak}.label];
        cfg = [];
        cfg.xlim = params.peaks{1}.peak_latency;
        %cfg.zlim = [0 6e-14];
        cfg.layout = 'neuromag306planar.lay';
        h = figure; ft_topoplotER(cfg,grandavg_squidgrad); colorbar; title(['GRAND AVG SQUID-GRAD - ' params.trigger_labels{i_trigger}])
        saveas(h, fullfile(base_save_path, 'figs', ['squidgrad' peak_label '_grndAvg_topo_trig-' params.trigger_labels{i_trigger} '.jpg'])) 
        close all
    end

    % SQUID-MAG
    cfg = [];
    cfg.channel = 'megmag';
    grandavg_squidmag = ft_timelockgrandaverage(cfg,squid{i_trigger,~cellfun(@isempty,squid(i_trigger,:))});

    h = figure; 
    plot(grandavg_squidmag.time*1e3,grandavg_squidmag.avg*1e15)
    xlabel('t [msec]')
    ylabel('B [fT]')
    xlim([-params.pre params.post]*1e3);
    title(['Grand average squidmag - phalange ' params.trigger_labels{i_trigger}])
    saveas(h, fullfile(base_save_path, 'figs', ['squidmag_grndAvg_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']))
    close all

    for i_peak = 1:length(params.peaks)
        peak_label = ['_' params.peaks{i_peak}.label];
        cfg = [];
        cfg.xlim = params.peaks{1}.peak_latency;
        %cfg.zlim = [0 6e-14];
        cfg.layout = 'neuromag306mag.lay';
        h = figure; ft_topoplotER(cfg,grandavg_squidmag); colorbar; title(['GRAND AVG SQUID-MAG - ' params.trigger_labels{i_trigger}])
        saveas(h, fullfile(base_save_path, 'figs', ['squidmag' peak_label '_grndAvg_topo_trig-' params.trigger_labels{i_trigger} '.jpg']))
        close all
    end

    clear grandavg_squidmag grandavg_squidgrad grandavg_opm
end
end
