function dipole_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_subs = max(subs);
n_ph = length(params.phalange_labels);
dist_sqmag_opm = nan(n_subs,n_ph);
dist_sqgrad_opm = nan(n_subs,n_ph);
dist_sqmag_sqgrad = nan(n_subs,n_ph);
spread_opm = nan(n_subs,n_ph);
spread_squidmag = nan(n_subs,n_ph);
spread_squidgrad = nan(n_subs,n_ph);
mom_squidmag = nan(n_subs,n_ph);
mom_squidgrad = nan(n_subs,n_ph);
mom_opm = nan(n_subs,n_ph);
for i_sub = subs
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    ft_hastoolbox('mne', 1);
    save_path = fullfile(base_save_path,params.sub);
    clear squidmag_dipole squidgrad_dipole opm_dipole
    dipole_squidmag{i_sub} = load(fullfile(save_path, 'dipoles')).squidmag_dipole;
    dipole_squidgrad{i_sub} = load(fullfile(save_path, 'dipoles')).squidgrad_dipole;
    dipole_opm{i_sub} = load(fullfile(save_path, 'dipoles')).opm_dipole;
  
    % Metrics: 
    % - distance between dipoles for same phalange different systems
    % - over phalanges: average distance from mean location within distance
    pos_squidmag = zeros(n_ph,3);
    pos_squidgrad = zeros(n_ph,3);
    pos_opm = zeros(n_ph,3);    

    for i_phalange = 1:n_ph
        pos_squidmag(i_phalange,:) = dipole_squidmag{i_sub}{i_phalange}.dip.pos;
        pos_squidgrad(i_phalange,:) = dipole_squidgrad{i_sub}{i_phalange}.dip.pos;
        pos_opm(i_phalange,:) = dipole_opm{i_sub}{i_phalange}.dip.pos;

        mom_squidmag(i_sub,i_phalange) = max(vecnorm(dipole_squidmag{i_sub}{i_phalange}.dip.mom,2,1));
        mom_squidgrad(i_sub,i_phalange) = max(vecnorm(dipole_squidgrad{i_sub}{i_phalange}.dip.mom,2,1));
        mom_opm(i_sub,i_phalange) = max(vecnorm(dipole_opm{i_sub}{i_phalange}.dip.mom,2,1));

        dist_sqmag_opm(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqgrad_opm(i_sub,i_phalange) = 1e1*norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqmag_sqgrad(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));
    end
    D = pdist2(pos_opm,pos_opm);
    goods = ~any(D<50&D>0,2);
    if sum(goods)>=3 % remove up to 2 outliers
        center = mean(pos_opm(goods,:),1);
    else
        center = mean(pos_opm,1);
    end
    spread_opm(i_sub,:) = 1e1*vecnorm(pos_opm-center,2,2)';

    D = pdist2(pos_squidmag,pos_squidmag);
    goods = ~any(D<50&D>0,2);
    if sum(goods)>=3 % remove up to 2 outliers
        center = mean(pos_squidmag(goods,:),1);
    else
        center = mean(pos_squidmag,1);
    end
    spread_squidmag(i_sub,:) = 1e1*vecnorm(pos_opm-center,2,2)'; % mean distance from center of phalanges
    
    D = pdist2(pos_squidgrad,pos_squidgrad);
    goods = ~any(D<50&D>0,2);
    if sum(goods)>=3 % remove up to 2 outliers 
        center = mean(pos_squidgrad(goods,:),1);
    else
        center = mean(pos_squidgrad,1);
    end
    spread_squidgrad(i_sub,:) = 1e1*vecnorm(pos_opm-center,2,2)'; % mean distance from center of phalanges
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
title(['Dist SQMAG to OPM (mean = ' num2str(mean(mean(dist_sqmag_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidmag_to_opm_dist.jpg'))

h = figure('DefaultAxesFontSize',16);
bar(1:length(params.phalange_labels),mean(dist_sqgrad_opm,1,'omitnan'));
hold on
er = errorbar(1:5,mean(dist_sqgrad_opm,1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-min(dist_sqgrad_opm,[],1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-max(dist_sqgrad_opm,[],1,'omitnan'));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
er.CapSize = 30;
hold off
title(['Dist SQGRAD to OPM (mean = ' num2str(mean(mean(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidgrad_to_opm_dist.jpg'))
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
title(['Dist SQMAG to SQGRAD (mean = ' num2str(mean(mean(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
ylabel('Distance [mm]')
xlabel('Phalange')
xticklabels(params.phalange_labels)
saveas(h, fullfile(base_save_path, 'figs', 'dipole_squidmag_to_squidgrad_dist.jpg'))
close

for i_ph = 1:5
    h = figure('DefaultAxesFontSize',16);
    plot(subs,dist_sqmag_opm(subs,i_ph),'+-');
    hold on
    plot(subs,dist_sqgrad_opm(subs,i_ph),'x-');
    plot(subs,dist_sqmag_sqgrad(subs,i_ph),'*-');
    hold off
    title([params.phalange_labels{i_ph} ' - dipole distances over subjects'])
    ylabel('Distance [mm]')
    xlabel('Subjects')
    legend(['SQMAG-OPM   '; 'SQGRAD-OPM  '; 'SQMAG-SQGRAD'])
    saveas(h, fullfile(base_save_path, 'figs', ['dipole_dist_vs_sub-' params.phalange_labels{i_ph} '.jpg']))
    close
end

%% Plot spread squidmag vs opm
data1 = spread_squidmag;
data2 = spread_opm;
data3 = spread_squidgrad;
mean1 = mean(data1,1,'omitnan');
mean2 = mean(data2,1,'omitnan');
mean3 = mean(data3,1,'omitnan');
min1 = min(data1,[],1,'omitnan');
min2 = min(data2,[],1,'omitnan');
min3 = min(data3,[],1,'omitnan');
max1 = max(data1,[],1,'omitnan');
max2 = max(data2,[],1,'omitnan');
max3 = max(data3,[],1,'omitnan');
err1 = [mean1-min1; max1-mean1];
err2 = [mean2-min2; max2-mean2];
err3 = [mean3-min3; max3-mean3];

h = figure('DefaultAxesFontSize',16);
h.Position(3) = round(h.Position(3)*1.3);
bar(1:length(params.phalange_labels),[mean1; mean2; mean3]','grouped');
hold on
for k=1:length(params.phalange_labels)
    errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
end
p_values = zeros(1, 5);
for i = 1:5
    [~, p_values(i,1)] = ttest(data1(:, i), data2(:, i));
    [~, p_values(i,2)] = ttest(data2(:, i), data3(:, i));
    [~, p_values(i,3)] = ttest(data1(:, i), data3(:, i));
    sigstar({[i-0.22, i]}, p_values(i,1));
    sigstar({[i, i+0.22]}, p_values(i,2));
    sigstar({[i-0.22, i+0.22]}, p_values(i,3));
end
hold off
title('Group level M60 dipole spread')
ylabel('Dipoles spread [mm]')
legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
saveas(h, fullfile(base_save_path, 'figs', 'dipole_spread.jpg'))
close

data = {spread_squidmag, spread_opm, spread_squidgrad};
triggerLabels = params.phalange_labels;
yLabelStr = 'Dipole spread [mm]';
titleStr = ['Group level ' params.peaks{1}.label ' dipole spread - SQMAG vs OPM vs SQGRAD'];
save_path = fullfile(base_save_path, 'figs', 'dipole_spread_sqmag_opm_sqgrad_box.jpg');
pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

% data = {spread_squidmag, spread_opm};
% triggerLabels = params.phalange_labels;
% yLabelStr = 'Dipole spread';
% titleStr = ['Group level ' params.peaks{1}.label ' dipole spread - SQMAG vs OPM'];
% save_path = fullfile(base_save_path, 'figs', 'dipole_spread_sqmag_opm_box.jpg');
% pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
% 
% data = {spread_squidgrad, spread_opm};
% triggerLabels = params.phalange_labels;
% yLabelStr = 'Dipole spread';
% titleStr = ['Group level ' params.peaks{1}.label ' dipole spread - SQGRAD vs OPM'];
% save_path = fullfile(base_save_path, 'figs', 'dipole_spread_sqgrad_opm_box.jpg');
% pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

%% Peak mom
data = {1e9*1e-4*mom_squidmag, 1e9*1e-4*mom_opm, 1e9*1e-4*mom_squidgrad};
triggerLabels = params.phalange_labels;
yLabelStr = 'Peak moment [nAm]';
titleStr = ['Group level ' params.peaks{1}.label ' peak dipole moment - SQMAG vs OPM vs SQGRAD'];
save_path = fullfile(base_save_path, 'figs', 'dipole_mom_sqmag_opm_sqgrad_box.jpg');
pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

data = {1e9*1e-4*mom_squidmag, 1e9*1e-4*mom_opm};
triggerLabels = params.phalange_labels;
yLabelStr = 'Peak moment [nAm]';
titleStr = ['Group level ' params.peaks{1}.label ' peak dipole moment - SQMAG vs OPM'];
save_path = fullfile(base_save_path, 'figs', 'dipole_mom_sqmag_opm_box.jpg');
pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

disp('done')
end