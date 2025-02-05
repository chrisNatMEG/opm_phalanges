function mne_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_subs = max(subs);
n_ph = length(params.phalange_labels);
dist_sqmag_opm = nan(n_subs,n_ph);
dist_sqgrad_opm = nan(n_subs,n_ph);
dist_sqmag_sqgrad = nan(n_subs,n_ph);
dist_sqeeg_opmeeg = nan(n_subs,n_ph);
fahm_opm = nan(n_subs,n_ph);
fahm_squidmag = nan(n_subs,n_ph);
fahm_squidgrad = nan(n_subs,n_ph);
fahm_squideeg = nan(n_subs,n_ph);
fahm_opmeeg = nan(n_subs,n_ph);
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    load(fullfile(save_path, 'mnes')); 
    mne_squidmag{i_sub} = squidmag_mne_M100;
    mne_squidgrad{i_sub} = squidgrad_mne_M100;
    mne_opm{i_sub} = opm_mne_M100;
    mne_squideeg{i_sub} = squideeg_mne_M100;
    mne_opmeeg{i_sub} = opmeeg_mne_M100;
  
    % Metrics: 
    % - distance between mnes for same phalange different systems
    % - over phalanges: average distance from mean location within distance
    pos_squidmag = zeros(n_ph,3);
    pos_squidgrad = zeros(n_ph,3);
    pos_opm = zeros(n_ph,3);
    pos_squideeg = zeros(n_ph,3);
    pos_opmeeg = zeros(n_ph,3);
    for i_phalange = 1:n_ph
        pos_squidmag(i_phalange,:) = mne_squidmag{i_sub}{i_phalange}.peakloc;
        pos_squidgrad(i_phalange,:) = mne_squidgrad{i_sub}{i_phalange}.peakloc;
        pos_opm(i_phalange,:) = mne_opm{i_sub}{i_phalange}.peakloc;

        dist_sqmag_opm(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqgrad_opm(i_sub,i_phalange) = 1e1*norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqmag_sqgrad(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));
        if ~isempty(mne_squideeg{i_sub})
            pos_squideeg(i_phalange,:) = 1e1*mne_squideeg{i_sub}{i_phalange}.peakloc;
            pos_opmeeg(i_phalange,:) = 1e1*mne_opmeeg{i_sub}{i_phalange}.peakloc;
            dist_sqeeg_opmeeg(i_sub,i_phalange) = 1e1*norm(pos_squideeg(i_phalange,:)-pos_opmeeg(i_phalange,:));
        end
    end
    fahm_opm(i_sub,:) = mne_opm{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
    fahm_squidmag(i_sub,:) = mne_squidmag{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
    fahm_squidgrad(i_sub,:) = mne_squidgrad{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
    if ~isempty(mne_squideeg{i_sub})
        fahm_squideeg(i_sub,:) = mne_squideeg{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
        fahm_opmeeg(i_sub,:) = mne_opmeeg{i_sub}{i_phalange}.fahm; % mean distance from center of phalanges
    end
end

%% Plot distances
h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqmag_opm,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqmag_opm,1,'omitnan'), mean(dist_sqmag_opm,1,'omitnan')-min(dist_sqmag_opm,[],1,'omitnan'), mean(dist_sqmag_opm,1,'omitnan')-max(dist_sqmag_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQUID-MAG to OPM (mean = ' num2str(mean(mean(dist_sqmag_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'mne_squidmag_to_opm_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqgrad_opm,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqgrad_opm,1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-min(dist_sqgrad_opm,[],1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-max(dist_sqgrad_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQUID-GRAD to OPM (mean = ' num2str(mean(mean(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'mne_squidgrad_to_opm_dist.jpg'))
close

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqmag_sqgrad,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqmag_sqgrad,1,'omitnan'), mean(dist_sqmag_sqgrad,1,'omitnan')-min(dist_sqmag_sqgrad,[],1,'omitnan'), mean(dist_sqmag_sqgrad,1,'omitnan')-max(dist_sqmag_sqgrad,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQUID-MAG to SQUID-GRAD (mean = ' num2str(mean(mean(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'mne_squidmag_to_squidgrad_dist.jpg'))
close

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqeeg_opmeeg,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqeeg_opmeeg,1,'omitnan'), mean(dist_sqeeg_opmeeg,1,'omitnan')-min(dist_sqeeg_opmeeg,[],1,'omitnan'), mean(dist_sqeeg_opmeeg,1,'omitnan')-max(dist_sqeeg_opmeeg,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQUID-EEG to OPM-EEG (mean = ' num2str(mean(mean(dist_sqeeg_opmeeg,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'mne_squideeg_to_opmeeg_dist.jpg'))
close all

%% Plot FAHM squidmag vs opm
data1 = fahm_squidmag;
data2 = fahm_opm;
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
title('MNE: Group level M100 FAHM')
ylabel('M100 FAHM [m^2]')
xlabel('Phalange')
legend({'squidmag','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'mne_fahm_squidmag_opm.jpg'))

%% Plot FAHM squidgrad vs opm
data1 = fahm_squidgrad;
data2 = fahm_opm;
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
title('MNE: Group level M100 FAHM')
ylabel('M100 FAHM [m^2]')
xlabel('Phalange')
legend({'squidgrad','opm'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'mne_fahm_squidgrad_opm.jpg'))

%% Plot FAHM squidgrad vs opm
data1 = fahm_squidgrad;
data2 = fahm_squidmag;
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
title('MNE: Group level M100 FAHM')
ylabel('M100 FAHM [m^2]')
xlabel('Phalange')
legend({'squidgrad','squidmag'});
xticklabels(params.phalange_labels)
saveas(h, fullfile(save_path, 'figs', 'mne_fahm_squidgrad_squidmag.jpg'))

end