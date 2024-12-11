function dipole_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

peak_ratio = [];
snr = [];
latency = [];
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    load(fullfile(save_path, 'dipoles')); 
    dipole_squidmag{i_sub} = megmag_dipole;
    dipole_squidgrad{i_sub} = megplanar_dipole;
    dipole_opm{i_sub} = opm_dipole;
    dipole_squideeg{i_sub} = eeg_dipole;
    dipole_opmeeg{i_sub} = opmeeg_dipole;

    n_ph = length(params.phalange_labels);
    % Metrics: 
    % - distance between dipoles for same phalange different systems
    % - over phalanges: average distance from mean location within distance
    pos_squidmag = zeros(n_ph,3);
    pos_squidgrad = zeros(n_ph,3);
    pos_opm = zeros(n_ph,3);
    pos_squideeg = zeros(n_ph,3);
    pos_opmeeg = zeros(n_ph,3);
    for i_phalange = 1:n_ph
        pos_squidmag(i_phalange,:) = dipole_squidmag{i_sub}.pos;
        pos_squidgrad(i_phalange,:) = dipole_squidgrad{i_sub}.pos;
        pos_opm(i_phalange,:) = dipole_opm{i_sub}.pos;
        pos_squideeg(i_phalange,:) = dipole_squideeg{i_sub}.pos;
        pos_opmeeg(i_phalange,:) = dipole_opmeeg{i_sub}.pos;

        dist_sqmag_opm.meg(i_sub,i_phalange) = norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqgrad_opm.meg(i_sub,i_phalange) = norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqmag_sqgrad.meg(i_sub,i_phalange) = norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));
        dist_sqeeg_opmeeg.meg(i_sub,i_phalange) = norm(pos_squideeg(i_phalange,:)-pos_opmeeg(i_phalange,:));
    end
    meandist_opm(i_sub) = mean(vecnorm(pos_opm-mean(pos_opm,1),2,2)); % mean distance from center of phalanges
end

%% Save
save(fullfile(save_path, 'group_sensor'),"peak_ratio""snr","latency","-v7.3");

%% Plot ratio
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.meg,1));
hold on
er = errorbar(1:5,mean(peak_ratio.meg,1), mean(peak_ratio.meg,1)-min(peak_ratio.meg,[],1), mean(peak_ratio.meg,1)-max(peak_ratio.meg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amplitude ratio (mean = ' num2str(mean(mean(peak_ratio.meg))) ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Peak_amplitude_ratios_meg.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(peak_ratio.eeg,1));
hold on
er = errorbar(1:5,mean(peak_ratio.eeg,1), mean(peak_ratio.eeg,1)-min(peak_ratio.eeg,[],1), mean(peak_ratio.eeg,1)-max(peak_ratio.eeg,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 peak amplitude ratio (mean = ' num2str(mean(mean(peak_ratio.eeg))) ')'])
ylabel('OPMEEG/SQUIDEEG')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Peak_amplitude_ratios_eeg.jpg'))

%% Plot SNR - error
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_error,1));
hold on
er = errorbar(1:5,mean(snr.ratio_error,1), mean(snr.ratio_error,1)-min(snr.ratio_error,[],1), mean(snr.ratio_error,1)-max(snr.ratio_error,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{stderror} ratio (mean = ' num2str(mean(mean(snr.ratio_error))) ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_ratios_error.jpg'))

%% Plot SNR - prestim
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(snr.ratio_prestim,1));
hold on
er = errorbar(1:5,mean(snr.ratio_prestim,1), mean(snr.ratio_prestim,1)-min(snr.ratio_prestim,[],1), mean(snr.ratio_prestim,1)-max(snr.ratio_prestim,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 SNR_{prestim} ratio (mean = ' num2str(mean(mean(snr.ratio_prestim))) ')'])
ylabel('OPM/SQUID')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'SNR_ratios_prestim.jpg'))

%% Plot peak latency
tmp = latency.opm-latency.meg;
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(tmp,1));
hold on
er = errorbar(1:5,mean(tmp,1), mean(tmp,1)-min(tmp,[],1), mean(tmp,1)-max(tmp,[],1));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['M100 latency diff (opm mean = ' num2str(mean(mean(latency.opm))) '; meg mean = ' num2str(mean(mean(latency.meg))) ')'])
ylabel('t_{OPM}-t_{SQUID}')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'Latency.jpg'))
end