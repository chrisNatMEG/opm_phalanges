function sensor_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

peak_ratio = [];
snr = [];
latency = [];
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    clear peak_opm peak_opmeeg
    clear peak_squidmag peak_squideeg
    peak_opm = load(fullfile(save_path, [params.sub '_opm_' params.peaks{1}.label '.mat'])).peak; 
    peak_opmeeg = load(fullfile(save_path, [params.sub '_opmeeg_' params.peaks{1}.label '.mat'])).peak; 
    peak_squidmag = load(fullfile(save_path, [params.sub '_squidmag_' params.peaks{1}.label '.mat'])).peak; 
    peak_squideeg = load(fullfile(save_path, [params.sub '_squideeg_' params.peaks{1}.label '.mat'])).peak;

    clear squidmag_timelocked opm_timelocked
    squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked'])).timelocked; 
    opm_timelocked = load(fullfile(save_path, [params.sub '_opm_timelocked'])).timelocked; 
    meg_chs = find(contains(squidmag_timelocked{1}.label,'MEG'));
    opm_chs = find(contains(opm_timelocked{1}.label,'bz'));

    for i_phalange = 1:length(params.phalange_labels)
        peak_ratio.meg(i_sub,i_phalange) = peak_opm{i_phalange}.peak_amplitude/peak_squidmag{i_phalange}.peak_amplitude;
        peak_ratio.eeg(i_sub,i_phalange) = peak_opmeeg{i_phalange}.peak_amplitude/peak_squideeg{i_phalange}.peak_amplitude;
        snr.error_opm(i_sub,i_phalange) = peak_opm{i_phalange}.peak_amplitude/peak_opm{i_phalange}.std_error;
        snr.error_squidmag(i_sub,i_phalange) = peak_squidmag{i_phalange}.peak_amplitude/peak_squidmag{i_phalange}.std_error;
        snr.error_squideeg(i_sub,i_phalange) = peak_squideeg{i_phalange}.peak_amplitude/peak_squideeg{i_phalange}.std_error;
        snr.error_opmeeg(i_sub,i_phalange) = peak_opmeeg{i_phalange}.peak_amplitude/peak_opmeeg{i_phalange}.std_error;
        snr.prestim_opm(i_sub,i_phalange) = peak_opm{i_phalange}.peak_amplitude/peak_opm{i_phalange}.prestim_std;
        snr.prestim_squidmag(i_sub,i_phalange) = peak_squidmag{i_phalange}.peak_amplitude/peak_squidmag{i_phalange}.prestim_std;
        snr.prestim_opmeeg(i_sub,i_phalange) = peak_opmeeg{i_phalange}.peak_amplitude/peak_opmeeg{i_phalange}.prestim_std;
        snr.prestim_squideeg(i_sub,i_phalange) = peak_squideeg{i_phalange}.peak_amplitude/peak_squideeg{i_phalange}.prestim_std;
        snr.ratio_error(i_sub,i_phalange) = snr.error_opm(i_sub,i_phalange)/snr.error_squidmag(i_sub,i_phalange);
        snr.ratio_prestim(i_sub,i_phalange) = snr.prestim_opm(i_sub,i_phalange)/snr.prestim_squidmag(i_sub,i_phalange);
        latency.opm(i_sub,i_phalange) = peak_opm{i_phalange}.peak_latency;
        latency.squidmag(i_sub,i_phalange) = peak_squidmag{i_phalange}.peak_latency;
        latency.opmeeg(i_sub,i_phalange) = peak_opmeeg{i_phalange}.peak_latency;
        latency.squideeg(i_sub,i_phalange) = peak_squideeg{i_phalange}.peak_latency;
        amp.opm(i_sub,i_phalange) = peak_opm{i_phalange}.peak_amplitude;
        amp.squidmag(i_sub,i_phalange) = peak_squidmag{i_phalange}.peak_amplitude;
        amp.opmeeg(i_sub,i_phalange) = peak_opmeeg{i_phalange}.peak_amplitude;
        amp.squideeg(i_sub,i_phalange) = peak_squideeg{i_phalange}.peak_amplitude;

        h = figure;
        subplot(2,1,1)
        plot(squidmag_timelocked{i_phalange}.time*1e3,squidmag_timelocked{i_phalange}.avg(meg_chs(1:3:end),:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        xlim([-params.pre params.post])
        title(['Evoked SQUID MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(squidmag_timelocked{i_phalange}.cfg.trials)) ')'])
        subplot(2,1,2)
        plot(opm_timelocked{i_phalange}.time*1e3,opm_timelocked{i_phalange}.avg(opm_chs,:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        xlim([-params.pre params.post])
        title(['Evoked OPM MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(opm_timelocked{i_phalange}.cfg.trials)) ')'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_squidopm_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
        close all
    end
end

%% Save
save(fullfile(base_save_path, 'group_sensor'),"peak_ratio","snr","latency","amp","-v7.3");

%% Plot SNR vs subs
for i_ph = 1:5
    h = figure('DefaultAxesFontSize',16);
    plot(subs,snr.error_opm(:,i_ph))
    hold on
    plot(subs,snr.error_squidmag(:,i_ph),':')
    hold off
    title("SNR_{error} over subjects (:=squidmag)")
    ylabel("SNR")
    xlabel("sub")
    saveas(h, fullfile(base_save_path, 'figs', ['Peak_SNR_vs_subs' params.phalange_labels{i_ph} '.jpg']))
end

%% Plot amp vs subs
for i_ph = 1:5
    h = figure('DefaultAxesFontSize',16);
    plot(subs,amp.opm(:,i_ph)*1e15)
    hold on
    plot(subs,amp.squidmag(:,i_ph)*1e15,':')
    hold off
    title("M60 amp over subjects (:=squidmag)")
    ylabel("fT")
    xlabel("sub")
    saveas(h, fullfile(base_save_path, 'figs', ['Peak_amp_vs_subs' params.phalange_labels{i_ph} '.jpg']))
end

%% Plot ratio
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.meg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(peak_ratio.meg,1,'omitnan'), mean(peak_ratio.meg,1,'omitnan')-min(peak_ratio.meg,[],1,'omitnan'), mean(peak_ratio.meg,1,'omitnan')-max(peak_ratio.meg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title([params.peaks{1}.label ' peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.meg,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'Peak_amplitude_ratios_meg.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.eeg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(peak_ratio.eeg,1,'omitnan'), mean(peak_ratio.eeg,1,'omitnan')-min(peak_ratio.eeg,[],1,'omitnan'), mean(peak_ratio.eeg,1,'omitnan')-max(peak_ratio.eeg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title([params.peaks{1}.label ' peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.eeg,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPMEEG/SQUIDEEG')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'Peak_amplitude_ratios_eeg.jpg'))

%% Plot SNR - error
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_error,1,'omitnan'));
hold on
er = errorbar(1:5,mean(snr.ratio_error,1,'omitnan'), mean(snr.ratio_error,1,'omitnan')-min(snr.ratio_error,[],1,'omitnan'), mean(snr.ratio_error,1,'omitnan')-max(snr.ratio_error,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title([params.peaks{1}.label ' SNR_{stderror} ratio (mean = ' num2str(mean(mean(snr.ratio_error,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_ratios_error.jpg'))

%% Plot SNR - prestim
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_prestim,1,'omitnan'));
hold on
er = errorbar(1:5,mean(snr.ratio_prestim,1,'omitnan'), mean(snr.ratio_prestim,1,'omitnan')-min(snr.ratio_prestim,[],1,'omitnan'), mean(snr.ratio_prestim,1,'omitnan')-max(snr.ratio_prestim,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title([params.peaks{1}.label ' SNR_{prestim} ratio (mean = ' num2str(mean(mean(snr.ratio_prestim,'omitnan'),'omitnan'),'%.2f') ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_ratios_prestim.jpg'))

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
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title(['Group level ' params.peaks{1}.label ' amplitude'])
ylabel('Peak amplitude [fT]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'squidmag','opm'});
saveas(h, fullfile(base_save_path, 'figs', 'Amplitude_meg.jpg'))

% EEG
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
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title(['Group level ' params.peaks{1}.label ' amplitude'])
ylabel('Peak amplitude [uV]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'squideeg','opmeeg'});
saveas(h, fullfile(base_save_path, 'figs', 'Amplitude_eeg.jpg'))

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
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title(['Group level ' params.peaks{1}.label ' latency'])
ylabel('Latency [ms]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'squidmag','opm'},'Location','southeast');
saveas(h, fullfile(base_save_path, 'figs', 'Latency.jpg'))

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
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title('Group level SNR_{stderror}')
ylabel('SNR')
xlabel('Phalange')
legend({'squidmag','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_error.jpg'))

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
bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
end
sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
hold off
title('Group level SNR_{prestim}')
ylabel('SNR')
xlabel('Phalange')
legend({'squidmag','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_prestim.jpg'))

close all

%% Grand average
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);

    squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked'])).timelocked; 
    opm_timelocked = load(fullfile(save_path, [params.sub '_opm_timelocked'])).timelocked; 

    for i_ph = 1:length(params.phalange_labels)
        squid{i_ph}{i_sub} = squid_timelocked{i_ph};
        opm{i_ph}{i_sub} = opm_timelocked{i_ph};
    end
    clear squid_timelocked opm_timelocked
end
for i_ph = 1:length(params.phalange_labels)
    % OPM
    cfg = [];
    cfg.channel = '*bz';
    grandavg_opm{i_ph} = ft_timelockgrandaverage(cfg,opm{i_ph}{:});

    h = figure; 
    plot(grandavg_opm{i_ph}.time*1e3,grandavg_opm{i_ph}.avg*1e15)
    xlabel('t [msec]')
    ylabel('B [fT]')
    xlim([-params.pre params.post]*1e3);
    title(['Grand average opm - phalange ' params.phalange_labels{i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_opm_grndAvg_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
    close all

    cfg = [];
    cfg.xlim = params.peaks{1}.peak_latency;
    %cfg.zlim = [0 6e-14];
    cfg.layout = 'fieldlinebeta2bz_helmet.mat';
    h = figure; ft_topoplotER(cfg,grandavg_opm{i_ph}); colorbar; title(['GRAND AVG OPM - ' params.phalange_labels{i_ph}])
    clsoe all

    % SQUID-GRAD
    cfg = [];
    cfg.channel = 'meggrad';
    grandavg_squidgrad{i_ph} = ft_timelockgrandaverage(cfg,squid{i_ph}{:});

    h = figure; 
    plot(grandavg_squidgrad{i_ph}.time*1e3,grandavg_squidgrad{i_ph}.avg*1e15)
    xlabel('t [msec]')
    ylabel('B [fT]')
    xlim([-params.pre params.post]*1e3);
    title(['Grand average squidgrad - phalange ' params.phalange_labels{i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_squidgrad_grndAvg_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
    close all

    cfg = [];
    cfg.xlim = params.peaks{1}.peak_latency;
    %cfg.zlim = [0 6e-14];
    cfg.layout = 'neuromag306planar.lay';
    h = figure; ft_topoplotER(cfg,grandavg_squidgrad{i_ph}); colorbar; title(['GRAND AVG SQUID-GRAD - ' params.phalange_labels{i_ph}])
    close all

    % SQUID-MAG
    cfg = [];
    cfg.channel = 'megmag';
    grandavg_squidmag{i_ph} = ft_timelockgrandaverage(cfg,squidmag{i_ph}{:});

    h = figure; 
    plot(grandavg_squidmag{i_ph}.time*1e3,grandavg_squidmag{i_ph}.avg*1e15)
    xlabel('t [msec]')
    ylabel('B [fT]')
    xlim([-params.pre params.post]*1e3);
    title(['Grand average squidmag - phalange ' params.phalange_labels{i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_squidmag_grndAvg_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
    close all

    cfg = [];
    cfg.xlim = params.peaks{1}.peak_latency;
    %cfg.zlim = [0 6e-14];
    cfg.layout = 'neuromag306mag.lay';
    h = figure; ft_topoplotER(cfg,grandavg_squidmag{i_ph}); colorbar; title(['GRAND AVG SQUID-MAG - ' params.phalange_labels{i_ph}])
    close all
end
end