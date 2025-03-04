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
    load(fullfile(save_path, [params.sub '_opm_M100'])); 
    M100_opm{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_opmeeg_M100'])); 
    M100_opmeeg{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_squidmag_M100'])); 
    M100_squidmag{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_squidgrad_M100'])); 
    M100_squidgrad{i_sub} = M100;
    load(fullfile(save_path, [params.sub '_squideeg_M100']));
    M100_squideeg{i_sub} = M100;
    
    load(fullfile(save_path, [params.sub '_squidmag_timelocked'])); 
    squidmag_timelocked = timelocked;
    load(fullfile(save_path, [params.sub '_opm_timelocked'])); 
    meg_chs = find(contains(squidmag_timelocked{1}.label,'MEG'));
    opm_timelocked = timelocked;
    opm_chs = find(contains(opm_timelocked{1}.label,'bz'));

    for i_phalange = 1:length(params.phalange_labels)
        peak_ratio.meg(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.peak_amplitude/M100_squidmag{i_sub}{i_phalange}.peak_amplitude;
        peak_ratio.eeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.peak_amplitude/M100_squideeg{i_sub}{i_phalange}.peak_amplitude;
        snr.error_opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.peak_amplitude/M100_opm{i_sub}{i_phalange}.std_error;
        snr.error_squid(i_sub,i_phalange) = M100_squidmag{i_sub}{i_phalange}.peak_amplitude/M100_squidmag{i_sub}{i_phalange}.std_error;
        snr.error_squideeg(i_sub,i_phalange) = M100_squideeg{i_sub}{i_phalange}.peak_amplitude/M100_squideeg{i_sub}{i_phalange}.std_error;
        snr.error_opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.peak_amplitude/M100_opmeeg{i_sub}{i_phalange}.std_error;
        snr.prestim_opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.peak_amplitude/M100_opm{i_sub}{i_phalange}.prestim_std;
        snr.prestim_squid(i_sub,i_phalange) = M100_squidmag{i_sub}{i_phalange}.peak_amplitude/M100_squidmag{i_sub}{i_phalange}.prestim_std;
        snr.prestim_opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.peak_amplitude/M100_opmeeg{i_sub}{i_phalange}.prestim_std;
        snr.prestim_squideeg(i_sub,i_phalange) = M100_squideeg{i_sub}{i_phalange}.peak_amplitude/M100_squideeg{i_sub}{i_phalange}.prestim_std;
        snr.ratio_error(i_sub,i_phalange) = snr.error_opm(i_sub,i_phalange)/snr.error_squid(i_sub,i_phalange);
        snr.ratio_prestim(i_sub,i_phalange) = snr.prestim_opm(i_sub,i_phalange)/snr.prestim_squid(i_sub,i_phalange);
        latency.opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.peak_latency;
        latency.squidmag(i_sub,i_phalange) = M100_squidmag{i_sub}{i_phalange}.peak_latency;
        latency.squidgrad(i_sub,i_phalange) = M100_squidgrad{i_sub}{i_phalange}.peak_latency;
        latency.opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.peak_latency;
        latency.squideeg(i_sub,i_phalange) = M100_squideeg{i_sub}{i_phalange}.peak_latency;
        amp.opm(i_sub,i_phalange) = M100_opm{i_sub}{i_phalange}.peak_amplitude;
        amp.squidmag(i_sub,i_phalange) = M100_squidmag{i_sub}{i_phalange}.peak_amplitude;
        amp.squidgrad(i_sub,i_phalange) = M100_squidgrad{i_sub}{i_phalange}.peak_amplitude;
        amp.opmeeg(i_sub,i_phalange) = M100_opmeeg{i_sub}{i_phalange}.peak_amplitude;
        amp.squideeg(i_sub,i_phalange) = M100_squideeg{i_sub}{i_phalange}.peak_amplitude;

        h = figure;
        subplot(2,1,1)
        plot(squidmag_timelocked{i_phalange}.time*1e3,squidmag_timelocked{i_phalange}.avg(meg_chs(1:3:end),:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked SQUID MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(squidmag_timelocked{i_phalange}.cfg.trials)) ')'])
        subplot(2,1,2)
        plot(opm_timelocked{i_phalange}.time*1e3,opm_timelocked{i_phalange}.avg(opm_chs,:)*1e15)
        xlabel('t [msec]')
        ylabel('B [fT]')
        title(['Evoked OPM MAG - phalange ' params.phalange_labels{i_phalange} ' (n_{trls}=' num2str(length(opm_timelocked{i_phalange}.cfg.trials)) ')'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_squidopm_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']))
    end
    close all
end

%% Save
save(fullfile(base_save_path, 'group_sensor'),"peak_ratio","snr","latency","amp","-v7.3");

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
title(['M100 peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.meg,'omitnan'),'omitnan'),'%.2f') ')'])
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
title(['M100 peak amp ratio (mean = ' num2str(mean(mean(peak_ratio.eeg,'omitnan'),'omitnan'),'%.2f') ')'])
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
title(['M100 SNR_{stderror} ratio (mean = ' num2str(mean(mean(snr.ratio_error,'omitnan'),'omitnan'),'%.2f') ')'])
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
title(['M100 SNR_{prestim} ratio (mean = ' num2str(mean(mean(snr.ratio_prestim,'omitnan'),'omitnan'),'%.2f') ')'])
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
title('Group level M100 amplitude')
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
title('Group level M100 amplitude')
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
title('Group level M100 latency')
ylabel('Latency [ms]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
legend({'squidmag','opm'},'Location','southeast');
saveas(h, fullfile(base_save_path, 'figs', 'Latency.jpg'))

%% Plot SNR - error
data1 = snr.error_squid;
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
title('Group level SNR_{m100,stderror}')
ylabel('SNR')
xlabel('Phalange')
legend({'squidmag','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_error.jpg'))

%% Plot SNR - prestim
data1 = snr.prestim_squid;
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
title('Group level SNR_{m100,prestim}')
ylabel('SNR')
xlabel('Phalange')
legend({'squidmag','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'SNR_prestim.jpg'))

close all
end