function dipole_results_goup(base_save_path, subs, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_subs = max(subs);

if params.numdipoles == 2
    dip_labels = {'_dip-L', '_dip-R'};
elseif params.numdipoles == 1
    dip_labels = {''};
end

for i_peak = 1:length(params.peaks)
    peak_label = ['_' params.peaks{i_peak}.label];
    for i_dip = 1:params.numdipoles
        n_triggers = length(params.trigger_labels);
        dist_sqmag_opm = nan(n_subs,n_triggers);
        dist_sqgrad_opm = nan(n_subs,n_triggers);
        dist_sqmag_sqgrad = nan(n_subs,n_triggers);
        distX_sqmag_opm = nan(n_subs,n_triggers);
        distX_sqgrad_opm = nan(n_subs,n_triggers);
        distX_sqmag_sqgrad = nan(n_subs,n_triggers);
        distY_sqmag_opm = nan(n_subs,n_triggers);
        distY_sqgrad_opm = nan(n_subs,n_triggers);
        distY_sqmag_sqgrad = nan(n_subs,n_triggers);
        distZ_sqmag_opm = nan(n_subs,n_triggers);
        distZ_sqgrad_opm = nan(n_subs,n_triggers);
        distZ_sqmag_sqgrad = nan(n_subs,n_triggers);
        spread_opm = nan(n_subs,n_triggers);
        spread_squidmag = nan(n_subs,n_triggers);
        spread_squidgrad = nan(n_subs,n_triggers);
        mom_squidmag = nan(n_subs,n_triggers);
        mom_squidgrad = nan(n_subs,n_triggers);
        mom_opm = nan(n_subs,n_triggers);
        dists_opm = nan(n_triggers,n_triggers,n_subs);
        dists_sqgrad = nan(n_triggers,n_triggers,n_subs);
        dists_sqmag = nan(n_triggers,n_triggers,n_subs);
        
        vec_sqmag_opm = nan(n_triggers,n_subs,3);
        vec_sqgrad_opm = nan(n_triggers,n_subs,3);
        vec_sqmag_sqgrad = nan(n_triggers,n_subs,3);
        vec_opm = nan(n_triggers,n_triggers,n_subs,3);
        vec_sqmag = nan(n_triggers,n_triggers,n_subs,3);
        vec_sqgrad = nan(n_triggers,n_triggers,n_subs,3);
        for i_sub = subs
            params.sub = ['sub_' num2str(i_sub,'%02d')];
            ft_hastoolbox('mne', 1);
            save_path = fullfile(base_save_path,params.sub);

            dipole_squidmag = load(fullfile(save_path,['squidmag_' params.peaks{i_peak}.label '_dipoles'])).dipoles;
            dipole_squidgrad = load(fullfile(save_path, ['squidgrad_' params.peaks{i_peak}.label '_dipoles'])).dipoles;
            dipole_opm = load(fullfile(save_path, ['opm_' params.peaks{i_peak}.label '_dipoles'])).dipoles;

            pos_squidmag = nan(n_triggers,3);
            pos_squidgrad = nan(n_triggers,3);
            pos_opm = nan(n_triggers,3);    
        
            for i_trigger = 1:n_triggers
                pos_squidmag(i_trigger,:) = dipole_squidmag{i_trigger}.dip.pos(i_dip,:);
                pos_squidgrad(i_trigger,:) = dipole_squidgrad{i_trigger}.dip.pos(i_dip,:);
                pos_opm(i_trigger,:) = dipole_opm{i_trigger}.dip.pos(i_dip,:);
        
                mom_squidmag(i_sub,i_trigger) = max(vecnorm(dipole_squidmag{i_trigger}.dip.mom(i_dip,:),2,1));
                mom_squidgrad(i_sub,i_trigger) = max(vecnorm(dipole_squidgrad{i_trigger}.dip.mom(i_dip,:),2,1));
                mom_opm(i_sub,i_trigger) = max(vecnorm(dipole_opm{i_trigger}.dip.mom(i_dip,:),2,1));
        
                dist_sqmag_opm(i_sub,i_trigger) = 1e1*norm(pos_squidmag(i_trigger,:)-pos_opm(i_trigger,:));
                dist_sqgrad_opm(i_sub,i_trigger) = 1e1*norm(pos_squidgrad(i_trigger,:)-pos_opm(i_trigger,:));
                dist_sqmag_sqgrad(i_sub,i_trigger) = 1e1*norm(pos_squidmag(i_trigger,:)-pos_squidgrad(i_trigger,:));
                distX_sqmag_opm(i_sub,i_trigger) = 1e1*(pos_squidmag(i_trigger,1)-pos_opm(i_trigger,1));
                distX_sqgrad_opm(i_sub,i_trigger) = 1e1*(pos_squidgrad(i_trigger,1)-pos_opm(i_trigger,1));
                distX_sqmag_sqgrad(i_sub,i_trigger) = 1e1*(pos_squidmag(i_trigger,1)-pos_squidgrad(i_trigger,1));
                distY_sqmag_opm(i_sub,i_trigger) = 1e1*(pos_squidmag(i_trigger,2)-pos_opm(i_trigger,2));
                distY_sqgrad_opm(i_sub,i_trigger) = 1e1*(pos_squidgrad(i_trigger,2)-pos_opm(i_trigger,2));
                distY_sqmag_sqgrad(i_sub,i_trigger) = 1e1*(pos_squidmag(i_trigger,2)-pos_squidgrad(i_trigger,2));
                distZ_sqmag_opm(i_sub,i_trigger) = 1e1*(pos_squidmag(i_trigger,3)-pos_opm(i_trigger,3));
                distZ_sqgrad_opm(i_sub,i_trigger) = 1e1*(pos_squidgrad(i_trigger,3)-pos_opm(i_trigger,3));
                distZ_sqmag_sqgrad(i_sub,i_trigger) = 1e1*(pos_squidmag(i_trigger,3)-pos_squidgrad(i_trigger,3));

                vec_sqmag_opm(i_trigger,i_sub,:) = 1e1*(pos_squidmag(i_trigger,:)-pos_opm(i_trigger,:));
                vec_sqgrad_opm(i_trigger,i_sub,:) = 1e1*(pos_squidgrad(i_trigger,:)-pos_opm(i_trigger,:));
                vec_sqmag_sqgrad(i_trigger,i_sub,:) = 1e1*(pos_squidmag(i_trigger,:)-pos_squidgrad(i_trigger,:));
                
                for j = 1:n_triggers
                    dists_opm(i_trigger,j,i_sub) = 1e1*norm(dipole_opm{i_trigger}.dip.pos(i_dip,:)-dipole_opm{j}.dip.pos(i_dip,:));
                    dists_sqmag(i_trigger,j,i_sub) = 1e1*norm(dipole_squidmag{i_trigger}.dip.pos(i_dip,:)-dipole_squidmag{j}.dip.pos(i_dip,:));
                    dists_sqgrad(i_trigger,j,i_sub) = 1e1*norm(dipole_squidgrad{i_trigger}.dip.pos(i_dip,:)-dipole_squidgrad{j}.dip.pos(i_dip,:));
                    vec_opm(i_trigger,j,i_sub,:) = 1e1*(dipole_opm{i_trigger}.dip.pos(i_dip,:)-dipole_opm{j}.dip.pos(i_dip,:));
                    vec_sqmag(i_trigger,j,i_sub,:) = 1e1*(dipole_squidmag{i_trigger}.dip.pos(i_dip,:)-dipole_squidmag{j}.dip.pos(i_dip,:));
                    vec_sqgrad(i_trigger,j,i_sub,:) = 1e1*(dipole_squidgrad{i_trigger}.dip.pos(i_dip,:)-dipole_squidgrad{j}.dip.pos(i_dip,:));
                end
            end
            clear squidmag_dipole squidgrad_dipole opm_dipole

            center = mean(pos_opm,1);
            spread_opm(i_sub,:) = 1e1*vecnorm(pos_opm-center,2,2)';
        
            center = mean(pos_squidmag,1);
            spread_squidmag(i_sub,:) = 1e1*vecnorm(pos_squidmag-center,2,2)'; % mean distance from center of phalanges
            
            center = mean(pos_squidgrad,1);
            spread_squidgrad(i_sub,:) = 1e1*vecnorm(pos_squidgrad-center,2,2)'; % mean distance from center of phalanges
        end

        %%
        save_path = fullfile(base_save_path, 'figs', ['dipole_dist_opm' peak_label dip_labels{i_dip} '.jpg']);
        plotDipDistances(dists_opm,'OPM',params,save_path,0)

        save_path = fullfile(base_save_path, 'figs', ['dipole_dist_sqmag' peak_label dip_labels{i_dip} '.jpg']);
        plotDipDistances(dists_sqmag,'SQMAG',params,save_path,0)

        save_path = fullfile(base_save_path, 'figs', ['dipole_dist_sqgrad' peak_label dip_labels{i_dip} '.jpg']);
        plotDipDistances(dists_sqgrad,'SQGRAD',params,save_path,0)
        
        save_path = fullfile(base_save_path, 'figs', ['dipole_vec_opm' peak_label dip_labels{i_dip} '.jpg']);
        plotDipDistanceVectors(vec_opm,'OPM',params,save_path,1)

        save_path = fullfile(base_save_path, 'figs', ['dipole_vec_sqmag' peak_label dip_labels{i_dip} '.jpg']);
        plotDipDistanceVectors(vec_sqmag,'SQMAG',params,save_path,1)

        save_path = fullfile(base_save_path, 'figs', ['dipole_vec_sqgrad' peak_label dip_labels{i_dip} '.jpg']);
        plotDipDistanceVectors(vec_sqgrad,'SQGRAD',params,save_path,1)

        %% Plot distances
        h = figure('DefaultAxesFontSize',16);
        bar(1:length(params.trigger_labels),median(dist_sqmag_opm,1,'omitnan'));
        hold on
        er = errorbar(1:n_triggers,median(dist_sqmag_opm,1,'omitnan'), median(dist_sqmag_opm,1,'omitnan')-min(dist_sqmag_opm,[],1,'omitnan'), median(dist_sqmag_opm,1,'omitnan')-max(dist_sqmag_opm,[],1,'omitnan'));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        er.LineWidth = 1;
        er.CapSize = 30;
        hold off
        title(['Dist SQMAG to OPM (mean = ' num2str(mean(median(dist_sqmag_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
        ylabel('Distance [mm]')
        xlabel('Phalange')
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['dipole_squidmag_to_opm_dist' peak_label dip_labels{i_dip} '.jpg']))
        
        h = figure('DefaultAxesFontSize',16);
        bar(1:length(params.trigger_labels),median(dist_sqgrad_opm,1,'omitnan'));
        hold on
        er = errorbar(1:n_triggers,median(dist_sqgrad_opm,1,'omitnan'), median(dist_sqgrad_opm,1,'omitnan')-min(dist_sqgrad_opm,[],1,'omitnan'), median(dist_sqgrad_opm,1,'omitnan')-max(dist_sqgrad_opm,[],1,'omitnan'));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        er.LineWidth = 1;
        er.CapSize = 30;
        hold off
        title(['Dist SQGRAD to OPM (mean = ' num2str(mean(median(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
        ylabel('Distance [mm]')
        xlabel('Phalange')
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['dipole_squidgrad_to_opm_dist' peak_label dip_labels{i_dip} '.jpg']))
        close
        
        h = figure('DefaultAxesFontSize',16);
        bar(1:length(params.trigger_labels),median(dist_sqmag_sqgrad,1,'omitnan'));
        hold on
        er = errorbar(1:n_triggers,median(dist_sqmag_sqgrad,1,'omitnan'), median(dist_sqmag_sqgrad,1,'omitnan')-min(dist_sqmag_sqgrad,[],1,'omitnan'), median(dist_sqmag_sqgrad,1,'omitnan')-max(dist_sqmag_sqgrad,[],1,'omitnan'));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        er.LineWidth = 1;
        er.CapSize = 30;
        hold off
        title(['Dist SQMAG to SQGRAD (mean = ' num2str(mean(median(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
        ylabel('Distance [mm]')
        xlabel('Phalange')
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['dipole_squidmag_to_squidgrad_dist' peak_label dip_labels{i_dip} '.jpg']))
        close
        
    %     for i_ph = 1:n_triggers
    %         h = figure('DefaultAxesFontSize',16);
    %         plot(subs,dist_sqmag_opm(subs,i_ph),'+-');
    %         hold on
    %         plot(subs,dist_sqgrad_opm(subs,i_ph),'x-');
    %         plot(subs,dist_sqmag_sqgrad(subs,i_ph),'*-');
    %         hold off
    %         title([params.trigger_labels{i_ph} ' - dipole distances over subjects'])
    %         ylabel('Distance [mm]')
    %         xlabel('Subjects')
    %         legend(['SQMAG-OPM   '; 'SQGRAD-OPM  '; 'SQMAG-SQGRAD'])
    %         saveas(h, fullfile(base_save_path, 'figs', ['dipole_dist_vs_sub-' params.trigger_labels{i_ph} dip_labels{i_dip} '.jpg']))
    %         close
    %     end
        
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
        bar(1:length(params.trigger_labels),[mean1; mean2; mean3]','grouped');
        hold on
        for k=1:length(params.trigger_labels)
            errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
            errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
            errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
        end
        p_values = zeros(1, 5);
        for i = 1:n_triggers
            [~, p_values(i,1)] = ttest(data1(:, i), data2(:, i));
            [~, p_values(i,2)] = ttest(data2(:, i), data3(:, i));
            [~, p_values(i,3)] = ttest(data1(:, i), data3(:, i));
            sigstar({[i-0.22, i]}, p_values(i,1));
            sigstar({[i, i+0.22]}, p_values(i,2));
            sigstar({[i-0.22, i+0.22]}, p_values(i,3));
        end
        hold off
        title(['Group level dipole spread (mean: ' num2str(mean(median(spread_squidmag,1,'omitnan')),'%.1f') ', ' num2str(mean(median(spread_opm,1,'omitnan')),'%.1f') ', ' num2str(mean(median(spread_squidgrad,1,'omitnan')),'%.1f') ')'])
        ylabel('Dipoles spread [mm]')
        legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
        saveas(h, fullfile(base_save_path, 'figs', ['dipole_spread' peak_label dip_labels{i_dip} '.jpg']))
        close
        
        data = {spread_squidmag, spread_opm, spread_squidgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Dipole spread [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' dipole spread - SQMAG vs OPM vs SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['dipole_spread_sqmag_opm_sqgrad' peak_label dip_labels{i_dip} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
        
        %% Dip distance
        data = {dist_sqmag_opm, dist_sqgrad_opm, dist_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Distance [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' dipole distance - SQMAG-OPM, SQGRAD-OPM, SQMAG-SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['dipole_dist' peak_label dip_labels{i_dip} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

        data = {distX_sqmag_opm, distX_sqgrad_opm, distX_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Distance [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' dipole distance X - SQMAG-OPM, SQGRAD-OPM, SQMAG-SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['dipole_distX' peak_label dip_labels{i_dip} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

        data = {distY_sqmag_opm, distY_sqgrad_opm, distY_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Distance [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' dipole distance Y - SQMAG-OPM, SQGRAD-OPM, SQMAG-SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['dipole_distY' peak_label dip_labels{i_dip} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

        data = {distZ_sqmag_opm, distZ_sqgrad_opm, distZ_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Distance [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' dipole distance Z - SQMAG-OPM, SQGRAD-OPM, SQMAG-SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['dipole_distZ' peak_label dip_labels{i_dip} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
  
        %% Peak mom
        data = {1e9*1e-4*mom_squidmag, 1e9*1e-4*mom_opm, 1e9*1e-4*1e-2*mom_squidgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak moment [nAm]';
        titleStr = ['Group level ' params.peaks{1}.label ' peak dipole moment - SQMAG vs OPM vs SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['dipole_mom_sqmag_opm_sqgrad' peak_label dip_labels{i_dip} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
        
        data = {1e9*1e-4*mom_squidmag, 1e9*1e-4*mom_opm};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak moment [nAm]';
        titleStr = ['Group level ' params.peaks{1}.label ' peak dipole moment - SQMAG vs OPM'];
        save_path = fullfile(base_save_path, 'figs', ['dipole_mom_sqmag_opm' peak_label dip_labels{i_dip} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

        %% Save
        mom = [];
        mom.opm = mom_opm;
        mom.sqmag = mom_squidmag;
        mom.sqgrad = mom_squidgrad;
        spread = [];
        spread.opm = spread_opm;
        spread.sqmag = spread_squidmag;
        spread.sqgrad = spread_squidgrad;
        save(fullfile(base_save_path, ['group_dip' peak_label dip_labels{i_dip}]),"dist_sqmag_opm","dist_sqgrad_opm","dist_sqmag_sqgrad","mom","spread","-v7.3");
    
        close all
    end
end
disp('done')
end