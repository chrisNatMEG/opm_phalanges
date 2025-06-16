function dipole_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_subs = max(subs);
n_ph = length(params.phalange_labels);
dist_sqmag_opm = nan(n_subs,n_ph);
dist_sqgrad_opm = nan(n_subs,n_ph);
dist_sqmag_sqgrad = nan(n_subs,n_ph);
spread_opm = nan(n_subs,1);
spread_squidmag = nan(n_subs,1);
spread_squidgrad = nan(n_subs,1);
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

        dist_sqmag_opm(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqgrad_opm(i_sub,i_phalange) = 1e1*norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
        dist_sqmag_sqgrad(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));
    end
    spread_opm(i_sub,:) = mean(1e1*vecnorm(pos_opm-repmat(mean(pos_opm,1),[n_ph 1]),2,2))'; % mean distance from center of phalanges
    spread_squidmag(i_sub,:) = mean(1e1*vecnorm(pos_squidmag-repmat(mean(pos_squidmag,1),[n_ph 1]),2,2))'; % mean distance from center of phalanges
    spread_squidgrad(i_sub,:) = mean(1e1*vecnorm(pos_squidgrad-repmat(mean(pos_squidgrad,1),[n_ph 1]),2,2))'; % mean distance from center of phalanges
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
title(['Dist SQUID-GRAD to OPM (mean = ' num2str(mean(mean(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
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
title(['Dist SQUID-MAG to SQUID-GRAD (mean = ' num2str(mean(mean(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
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
bar(1,[mean1; mean2; mean3]','grouped');
hold on
errorbar(1-0.22,mean1(1),err1(1,1),err1(2,1),'k','linestyle','none');
errorbar(1,mean2(1),err2(1,1),err2(2,1),'k','linestyle','none');
errorbar(1+0.22,mean3(1),err3(1,1),err3(2,1),'k','linestyle','none');
p_values = zeros(1, 3);
[~, p_values(1)] = ttest(data1, data2);
[~, p_values(2)] = ttest(data2, data3);
[~, p_values(3)] = ttest(data1, data3);
sigstar({[1-0.22, 1]}, p_values(1));
sigstar({[1, 1+0.22]}, p_values(2));
sigstar({[1-0.22, 1+0.22]}, p_values(3));
hold off
title('Group level M60 dipole spread')
ylabel('Dipoles spread [mm]')
legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
saveas(h, fullfile(base_save_path, 'figs', 'dipole_spread.jpg'))

close all
end