
function mne_results_goup(base_save_path, subs, sourcemodel, sourcemodel_inflated, params)
n_triggers = length(params.trigger_codes);
n_subs = length(subs);
cov = '';
if isfield(params,'use_cov') && strcmp(params.use_cov,'all')
    cov = '_covAll';
elseif isfield(params,'use_cov') && strcmp(params.use_cov,'resting_state')
    cov = '_covRS';
elseif isfield(params,'use_cov') && strcmp(params.use_cov,'empty_room')
    cov = '_covER';
end

if params.numdipoles == 2
    hemi_labels = {'_RH', '_LH'};
elseif params.numdipoles == 1
    hemi_labels = {''};
end

for i_peak = 1:length(params.peaks)
    peak_label = ['_' params.peaks{i_peak}.label];
    for i_hemi = 1:params.numdipoles
        i_sub = subs(1);
        params.sub = ['sub_' num2str(i_sub,'%02d')];
        ft_hastoolbox('mne', 1);
        save_path = fullfile(base_save_path,params.sub);
        tmp = load(fullfile(save_path,'mne_distributions.mat'));
        squidmag_peak = load(fullfile(save_path,'squidmag_mne_peaks.mat')).peaks;
        squidgrad_peak = load(fullfile(save_path,'squidgrad_mne_peaks.mat')).peaks;
        opm_peak = load(fullfile(save_path,'opm_mne_peaks.mat')).peaks;
        for i_trigger = 1:n_triggers
            squidmag_mne{i_trigger} = tmp.squidmag_mne{i_trigger};
            squidgrad_mne{i_trigger} = tmp.squidgrad_mne{i_trigger};
            opm_mne{i_trigger} = tmp.opm_mne{i_trigger};
            squidmag_mne{i_trigger}.avg.pow = tmp.squidmag_mne{i_trigger}.avg.pow/n_subs;
            squidgrad_mne{i_trigger}.avg.pow = tmp.squidgrad_mne{i_trigger}.avg.pow/n_subs;
            opm_mne{i_trigger}.avg.pow = tmp.opm_mne{i_trigger}.avg.pow/n_subs;
            squidmag_peak{i_trigger}.fahm = squidmag_peak{i_trigger}.fahm/n_subs;
            squidgrad_peak{i_trigger}.fahm = squidgrad_peak{i_trigger}.fahm/n_subs;
            opm_peak{i_trigger}.fahm = opm_peak{i_trigger}.fahm/n_subs;
            squidmag_peak{i_trigger}.latency = squidmag_peak{i_trigger}.latency/n_subs;
            squidgrad_peak{i_trigger}.latency = squidgrad_peak{i_trigger}.latency/n_subs;
            opm_peak{i_trigger}.latency = opm_peak{i_trigger}.latency/n_subs;
        end
        clear tmp
        for i_sub = subs(2:end)
            params.sub = ['sub_' num2str(i_sub,'%02d')];
            ft_hastoolbox('mne', 1);
            save_path = fullfile(base_save_path,params.sub);
            tmp = load(fullfile(save_path,'mne_distributions.mat'));
            t_squidmag_peak = load(fullfile(save_path,'squidmag_mne_peaks.mat')).peaks;
            t_squidgrad_peak = load(fullfile(save_path,'squidgrad_mne_peaks.mat')).peaks;
            t_opm_peak = load(fullfile(save_path,'opm_mne_peaks.mat')).peaks;
            for i_trigger = 1:n_triggers
                squidmag_mne{i_trigger}.avg.pow = squidmag_mne{i_trigger}.avg.pow + tmp.squidmag_mne{i_trigger}.avg.pow/n_subs;
                squidgrad_mne{i_trigger}.avg.pow = squidgrad_mne{i_trigger}.avg.pow + tmp.squidgrad_mne{i_trigger}.avg.pow/n_subs;
                opm_mne{i_trigger}.avg.pow = opm_mne{i_trigger}.avg.pow + tmp.opm_mne{i_trigger}.avg.pow/n_subs;
                squidmag_peak{i_trigger}.fahm = squidmag_peak{i_trigger}.fahm + t_squidmag_peak{i_trigger}.fahm/n_subs;
                squidgrad_peak{i_trigger}.fahm = squidgrad_peak{i_trigger}.fahm + t_squidgrad_peak{i_trigger}.fahm/n_subs;
                opm_peak{i_trigger}.fahm = opm_peak{i_trigger}.fahm + t_opm_peak{i_trigger}.fahm/n_subs;
                squidmag_peak{i_trigger}.latency = squidmag_peak{i_trigger}.latency + t_squidmag_peak{i_trigger}.latency/n_subs;
                squidgrad_peak{i_trigger}.latency = squidgrad_peak{i_trigger}.latency + t_squidgrad_peak{i_trigger}.latency/n_subs;
                opm_peak{i_trigger}.latency = opm_peak{i_trigger}.latency + t_opm_peak{i_trigger}.latency/n_subs;
            end
            clear tmp
        end
        for i_trigger = 1:n_triggers
            opm_mne{i_trigger}.pos = sourcemodel_inflated.pos;
            opm_mne{i_trigger}.tri = sourcemodel_inflated.tri;
            params.modality = 'opm';
            h =  plot_source_distribution(opm_mne{i_trigger},opm_peak{i_trigger},params,0);
            saveas(h, fullfile(base_save_path, 'figs', ['mne_grnd_avg_infl_opm_' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
            close all
            opm_mne{i_trigger}.pos = sourcemodel.pos;
            opm_mne{i_trigger}.tri = sourcemodel.tri;
            params.modality = 'opm';
            h =  plot_source_distribution(opm_mne{i_trigger},opm_peak{i_trigger},params,0);
            saveas(h, fullfile(base_save_path, 'figs', ['mne_grnd_avg_opm_' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
            close all

            squidmag_mne{i_trigger}.pos = sourcemodel_inflated.pos;
            squidmag_mne{i_trigger}.tri = sourcemodel_inflated.tri;
            params.modality = 'sqmag';
            h =  plot_source_distribution(squidmag_mne{i_trigger},squidmag_peak{i_trigger},params,0);
            saveas(h, fullfile(base_save_path, 'figs', ['mne_grnd_avg_infl_squidmag_' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
            close all
            squidmag_mne{i_trigger}.pos = sourcemodel.pos;
            squidmag_mne{i_trigger}.tri = sourcemodel.tri;
            params.modality = 'sqmag';
            h =  plot_source_distribution(squidmag_mne{i_trigger},squidmag_peak{i_trigger},params,0);
            saveas(h, fullfile(base_save_path, 'figs', ['mne_grnd_avg_squidmag_' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
            close all

            squidgrad_mne{i_trigger}.pos = sourcemodel_inflated.pos;
            squidgrad_mne{i_trigger}.tri = sourcemodel_inflated.tri;
            params.modality = 'sqgrad';
            h =  plot_source_distribution(squidgrad_mne{i_trigger},squidgrad_peak{i_trigger},params,0);
            saveas(h, fullfile(base_save_path, 'figs', ['mne_grnd_avg_infl_squidgrad_' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
            close all
            squidgrad_mne{i_trigger}.pos = sourcemodel.pos;
            squidgrad_mne{i_trigger}.tri = sourcemodel.tri;
            params.modality = 'sqgrad';
            h =  plot_source_distribution(squidgrad_mne{i_trigger},squidgrad_peak{i_trigger},params,0);
            saveas(h, fullfile(base_save_path, 'figs', ['mne_grnd_avg_squidgrad_' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
            close all
        end
    end
end
%%
for i_peak = 1:length(params.peaks)
    peak_label = ['_' params.peaks{i_peak}.label];
    for i_hemi = 1:params.numdipoles
        n_subs = max(subs);
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
        fahm_opm = nan(n_subs,n_triggers);
        fahm_squidmag = nan(n_subs,n_triggers);
        fahm_squidgrad = nan(n_subs,n_triggers);
        targetregion_opm = nan(n_subs,n_triggers);
        targetregion_squidmag = nan(n_subs,n_triggers);
        targetregion_squidgrad = nan(n_subs,n_triggers);
        lat_opm = nan(n_subs,n_triggers);
        lat_squidmag = nan(n_subs,n_triggers);
        lat_squidgrad = nan(n_subs,n_triggers);
        pow_opm = nan(n_subs,n_triggers);
        pow_squidmag = nan(n_subs,n_triggers);
        pow_squidgrad = nan(n_subs,n_triggers);
        mom_opm = nan(n_subs,n_triggers);
        mom_squidmag = nan(n_subs,n_triggers);
        mom_squidgrad = nan(n_subs,n_triggers);
        dists_opm = nan(n_triggers,n_triggers,n_subs);
        dists_sqgrad = nan(n_triggers,n_triggers,n_subs);
        dists_sqmag = nan(n_triggers,n_triggers,n_subs);
        vec_opm = nan(n_triggers,n_triggers,n_subs,3);
        vec_sqmag = nan(n_triggers,n_triggers,n_subs,3);
        vec_sqgrad = nan(n_triggers,n_triggers,n_subs,3);
        for i_sub = subs
            params.sub = ['sub_' num2str(i_sub,'%02d')];
            ft_hastoolbox('mne', 1);
            save_path = fullfile(base_save_path,params.sub);
            mne_squidmag = load(fullfile(save_path, 'squidmag_mne_peaks.mat')).peaks; 
            mne_squidgrad = load(fullfile(save_path, 'squidgrad_mne_peaks.mat')).peaks;
            mne_opm = load(fullfile(save_path, 'opm_mne_peaks.mat')).peaks;
          
            % Metrics: 
            % - distance between mnes for same phalange different systems
            % - over phalanges: average distance from mean location within distance
        
            for i_trigger = 1:n_triggers
                pos_squidmag = mne_squidmag{i_trigger,i_peak}.loc(i_hemi,:);
                pos_squidgrad = mne_squidgrad{i_trigger,i_peak}.loc(i_hemi,:);
                pos_opm = mne_opm{i_trigger,i_peak}.loc(i_hemi,:);
        
                dist_sqmag_opm(i_sub,i_trigger) = 1e1*norm(pos_squidmag-pos_opm);
                dist_sqgrad_opm(i_sub,i_trigger) = 1e1*norm(pos_squidgrad-pos_opm);
                dist_sqmag_sqgrad(i_sub,i_trigger) = 1e1*norm(pos_squidmag-pos_squidgrad);
                distX_sqmag_opm(i_sub,i_trigger) = 1e1*(pos_squidmag(1)-pos_opm(1));
                distX_sqgrad_opm(i_sub,i_trigger) = 1e1*(pos_squidgrad(1)-pos_opm(1));
                distX_sqmag_sqgrad(i_sub,i_trigger) = 1e1*(pos_squidmag(1)-pos_squidgrad(1));
                distY_sqmag_opm(i_sub,i_trigger) = 1e1*(pos_squidmag(2)-pos_opm(2));
                distY_sqgrad_opm(i_sub,i_trigger) = 1e1*(pos_squidgrad(2)-pos_opm(2));
                distY_sqmag_sqgrad(i_sub,i_trigger) = 1e1*(pos_squidmag(2)-pos_squidgrad(2));
                distZ_sqmag_opm(i_sub,i_trigger) = 1e1*(pos_squidmag(3)-pos_opm(3));
                distZ_sqgrad_opm(i_sub,i_trigger) = 1e1*(pos_squidgrad(3)-pos_opm(3));
                distZ_sqmag_sqgrad(i_sub,i_trigger) = 1e1*(pos_squidmag(3)-pos_squidgrad(3));

                for j = 1:n_triggers
                    dists_opm(i_trigger,j,i_sub) = 1e1*norm(mne_opm{i_trigger,i_peak}.loc(i_hemi,:)-mne_opm{j,i_peak}.loc(i_hemi,:));
                    dists_sqmag(i_trigger,j,i_sub) = 1e1*norm(mne_squidmag{i_trigger,i_peak}.loc(i_hemi,:)-mne_squidmag{j,i_peak}.loc(i_hemi,:));
                    dists_sqgrad(i_trigger,j,i_sub) = 1e1*norm(mne_squidgrad{i_trigger,i_peak}.loc(i_hemi,:)-mne_squidgrad{j,i_peak}.loc(i_hemi,:));
                    vec_opm(i_trigger,j,i_sub,:) = 1e1*(mne_opm{i_trigger,i_peak}.loc(i_hemi,:)-mne_opm{j,i_peak}.loc(i_hemi,:));
                    vec_sqmag(i_trigger,j,i_sub,:) = 1e1*(mne_squidmag{i_trigger,i_peak}.loc(i_hemi,:)-mne_squidmag{j,i_peak}.loc(i_hemi,:));
                    vec_sqgrad(i_trigger,j,i_sub,:) = 1e1*(mne_squidgrad{i_trigger,i_peak}.loc(i_hemi,:)-mne_squidgrad{j,i_peak}.loc(i_hemi,:));
                end
        
                pow_squidmag(i_sub,i_trigger) = mne_squidmag{i_trigger,i_peak}.pow(i_hemi);
                pow_squidgrad(i_sub,i_trigger) = mne_squidgrad{i_trigger,i_peak}.pow(i_hemi);
                pow_opm(i_sub,i_trigger) = mne_opm{i_trigger,i_peak}.pow(i_hemi);
        
                mom_squidmag(i_sub,i_trigger) = mne_squidmag{i_trigger,i_peak}.mom(i_hemi);
                mom_squidgrad(i_sub,i_trigger) = mne_squidgrad{i_trigger,i_peak}.mom(i_hemi);
                mom_opm(i_sub,i_trigger) = mne_opm{i_trigger,i_peak}.mom(i_hemi);
        
                lat_squidmag(i_sub,i_trigger) = mne_squidmag{i_trigger,i_peak}.latency;
                lat_squidgrad(i_sub,i_trigger) = mne_squidgrad{i_trigger,i_peak}.latency;
                lat_opm(i_sub,i_trigger) = mne_opm{i_trigger,i_peak}.latency;
        
                fahm_opm(i_sub,i_trigger) = mne_opm{i_trigger,i_peak}.fahm(i_hemi); % mean distance from center of phalanges
                fahm_squidmag(i_sub,i_trigger) = mne_squidmag{i_trigger,i_peak}.fahm(i_hemi); % mean distance from center of phalanges
                fahm_squidgrad(i_sub,i_trigger) = mne_squidgrad{i_trigger,i_peak}.fahm(i_hemi); % mean distance from center of phalanges
    
                targetregion_opm(i_sub,i_trigger) = mne_opm{i_trigger,i_peak}.target_region(i_hemi); % mean distance from center of phalanges
                targetregion_squidmag(i_sub,i_trigger) = mne_squidmag{i_trigger,i_peak}.target_region(i_hemi); % mean distance from center of phalanges
                targetregion_squidgrad(i_sub,i_trigger) = mne_squidgrad{i_trigger,i_peak}.target_region(i_hemi); % mean distance from center of phalanges
        
                overlap_opm_squidmag(i_sub,i_trigger) = mne_opm{i_trigger,i_peak}.overlap_squidmag(i_hemi)/mne_opm{i_trigger,i_peak}.fahm(i_hemi);
                overlap_opm_squidgrad(i_sub,i_trigger) = mne_opm{i_trigger,i_peak}.overlap_squidgrad(i_hemi)/mne_opm{i_trigger,i_peak}.fahm(i_hemi);
                overlap_squidmag_squidgrad(i_sub,i_trigger) = mne_squidmag{i_trigger,i_peak}.overlap_squidgrad(i_hemi)/mne_squidmag{i_trigger,i_peak}.fahm(i_hemi);
                overlap_squidmag_opm(i_sub,i_trigger) = mne_squidmag{i_trigger,i_peak}.overlap_opm(i_hemi)/mne_squidmag{i_trigger,i_peak}.fahm(i_hemi);
                overlap_squidgrad_squidmag(i_sub,i_trigger) = mne_squidgrad{i_trigger,i_peak}.overlap_squidmag(i_hemi)/mne_squidgrad{i_trigger,i_peak}.fahm(i_hemi);
                overlap_squidgrad_opm(i_sub,i_trigger) = mne_squidgrad{i_trigger,i_peak}.overlap_opm(i_hemi)/mne_squidgrad{i_trigger,i_peak}.fahm(i_hemi);
            end
            clear mne_opm mne_squidmag mne_squidgrad
        end

        %% Save
        target = [];
        target.opm = targetregion_opm;
        target.sqmag = targetregion_squidmag;
        target.sqgrad = targetregion_squidgrad;
        fahm = [];
        fahm.opm = fahm_opm;
        fahm.sqmag = fahm_squidmag;
        fahm.sqgrad = fahm_squidgrad;
        lat = [];
        lat.opm = lat_opm;
        lat.sqmag = lat_squidmag;
        lat.sqgrad = lat_squidgrad;
        save(fullfile(base_save_path, ['group_mne' peak_label hemi_labels{i_hemi}]),"target","fahm","lat","dist_sqmag_opm","dist_sqgrad_opm","dist_sqmag_sqgrad","-v7.3");

        %%
        save_path = fullfile(base_save_path, 'figs', ['mne_dist_opm' peak_label hemi_labels{i_hemi} '.jpg']);
        plotDipDistances(dists_opm,'OPM',params,save_path,0)

        save_path = fullfile(base_save_path, 'figs', ['mne_dist_sqmag' peak_label hemi_labels{i_hemi} '.jpg']);
        plotDipDistances(dists_sqmag,'SQMAG',params,save_path,0)

        save_path = fullfile(base_save_path, 'figs', ['mne_dist_sqgrad' peak_label hemi_labels{i_hemi} '.jpg']);
        plotDipDistances(dists_sqgrad,'SQGRAD',params,save_path,0)

        save_path = fullfile(base_save_path, 'figs', ['mne_vec_opm' peak_label hemi_labels{i_hemi} '.jpg']);
        plotDipDistanceVectors(vec_opm,'OPM',params,save_path,1)

        save_path = fullfile(base_save_path, 'figs', ['mne_vec_sqmag' peak_label hemi_labels{i_hemi} '.jpg']);
        plotDipDistanceVectors(vec_sqmag,'SQMAG',params,save_path,1)

        save_path = fullfile(base_save_path, 'figs', ['mne_vec_sqgrad' peak_label hemi_labels{i_hemi} '.jpg']);
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
        title(['Dist SQ-MAG to OPM (mean = ' num2str(mean(median(dist_sqmag_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
        ylabel('Distance [mm]')
        xlabel('Phalange')
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_squidmag_to_opm_dist_' peak_label cov hemi_labels{i_hemi} '.jpg']))
        
        h = figure('DefaultAxesFontSize',16);
        bar(1:length(params.trigger_labels),median(dist_sqgrad_opm,1,'omitnan'));
        hold on
        er = errorbar(1:n_triggers,median(dist_sqgrad_opm,1,'omitnan'), median(dist_sqgrad_opm,1,'omitnan')-min(dist_sqgrad_opm,[],1,'omitnan'), median(dist_sqgrad_opm,1,'omitnan')-max(dist_sqgrad_opm,[],1,'omitnan'));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        er.LineWidth = 1;
        er.CapSize = 30;
        hold off
        title(['Dist SQ-GRAD to OPM (mean = ' num2str(mean(median(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
        ylabel('Distance [mm]')
        xlabel('Phalange')
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_squidgrad_to_opm_dist_' peak_label cov hemi_labels{i_hemi} '.jpg']))
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
        title(['Dist SQ-MAG to SQ-GRAD (mean = ' num2str(mean(median(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
        ylabel('Distance [mm]')
        xlabel('Phalange')
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_squidmag_to_squidgrad_dist_' peak_label cov hemi_labels{i_hemi} '.jpg']))
        close
        
%         for i_trigger = 1:n_triggers
%             h = figure('DefaultAxesFontSize',16);
%             plot(subs,dist_sqmag_opm(subs,i_trigger));
%             hold on
%             plot(subs,dist_sqgrad_opm(subs,i_trigger));
%             plot(subs,dist_sqmag_sqgrad(subs,i_trigger));
%             hold off
%             title([params.trigger_labels{i_trigger} ' - mne distances over subjects'])
%             ylabel('Distance [mm]')
%             xlabel('Subjects')
%             legend(['SQMAG-OPM   '; 'SQGRAD-OPM  '; 'SQMAG-SQGRAD'])
%             saveas(h, fullfile(base_save_path, 'figs', ['mne_dist_vs_sub-' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
%             close
%         end
        
        %% Plot FAHM
        data1 = fahm_squidmag;
        data2 = fahm_opm;
        data3 = fahm_squidgrad;
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
        h.Position(3) = round(h.Position(3)*1.2);
        bar(1:length(params.trigger_labels),[mean1; mean2; mean3]','grouped');
        hold on
        for k=1:length(params.trigger_labels)
            errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
            errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
            errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
        end
        
        p_values = zeros(n_triggers, 3);
        for i = 1:n_triggers
            [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
            [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
            [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
        end
        for i = 1:n_triggers
            sigstar({[i-0.22, i]}, p_values(i, 1));
            sigstar({[i, i+0.22]}, p_values(i, 2));
            sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
        end
        
        hold off
        title(['MNE: Group level FAHM (mean: ' num2str(mean(median(fahm_squidmag,1,'omitnan')),'%.1f') ', ' num2str(mean(median(fahm_opm,1,'omitnan')),'%.1f') ', ' num2str(mean(median(fahm_squidgrad,1,'omitnan')),'%.1f') ')'])
        ylabel('FAHM [cm^2]')
        xlabel('Phalange')
        legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_fahm' peak_label cov hemi_labels{i_hemi} '.jpg']))
        close all

        data = {fahm_squidmag, fahm_opm, fahm_squidgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'FAHM [cm^2]';
        titleStr = ['Group level ' params.peaks{1}.label ' FAHM - SQM vs OPM vs SQG'];
        save_path = fullfile(base_save_path, 'figs', ['mne_fahm_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
        %% Plot peak powers
        data1 = pow_squidmag;
        data2 = pow_opm;
        data3 = pow_squidgrad;
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
        h.Position(3) = round(h.Position(3)*1.2);
        bar(1:length(params.trigger_labels),[mean1; mean2; mean3]','grouped');
        hold on
        for k=1:length(params.trigger_labels)
            errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
            errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
            errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
        end
        
        p_values = zeros(n_triggers, 3);
        for i = 1:n_triggers
            [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
            [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
            [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
        end
        for i = 1:n_triggers
            sigstar({[i-0.22, i]}, p_values(i, 1));
            sigstar({[i, i+0.22]}, p_values(i, 2));
            sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
        end
        
        hold off
        title('MNE: Group level peak power')
        ylabel('M60 pow')
        xlabel('Phalange')
        legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_pow_' peak_label cov hemi_labels{i_hemi} '.jpg']))
        
        data = {pow_squidmag, pow_opm, pow_squidgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak power';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Power - SQMAG vs OPM vs SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['mne_pow_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
        %% Plot peak mom
        data1 = mom_squidmag;
        data2 = mom_opm;
        data3 = mom_squidgrad;
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
        h.Position(3) = round(h.Position(3)*1.2);
        bar(1:length(params.trigger_labels),[mean1; mean2; mean3]','grouped');
        hold on
        for k=1:length(params.trigger_labels)
            errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
            errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
            errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
        end
        
        p_values = zeros(n_triggers, 3);
        for i = 1:n_triggers
            [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
            [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
            [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
        end
        for i = 1:n_triggers
            sigstar({[i-0.22, i]}, p_values(i, 1));
            sigstar({[i, i+0.22]}, p_values(i, 2));
            sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
        end
        
        hold off
        title('MNE: Group level peak moment')
        ylabel('M60 mom')
        xlabel('Phalange')
        legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_mom_' peak_label cov hemi_labels{i_hemi} '.jpg']))
        
        data = {1e9*1e-4*mom_squidmag, 1e9*1e-4*mom_opm, 1e9*1e-4*1e-2*mom_squidgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak moment [nAm]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Moment - SQMAG vs OPM vs SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['mne_mom_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
        
        data = {1e9*1e-4*mom_squidmag, 1e9*1e-4*mom_opm};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak moment [nAm]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Moment - SQMAG vs OPM'];
        save_path = fullfile(base_save_path, 'figs', ['mne_mom_sqmag_opm' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
        %% Plot peak powers - opm and squidmag only
        data1 = pow_squidmag;
        data2 = pow_opm;
        mean1 = mean(data1,1,'omitnan');
        mean2 = mean(data2,1,'omitnan');
        min1 = min(data1,[],1,'omitnan');
        min2 = min(data2,[],1,'omitnan');
        max1 = max(data1,[],1,'omitnan');
        max2 = max(data2,[],1,'omitnan');
        err1 = [mean1-min1; max1-mean1];
        err2 = [mean2-min2; max2-mean2];
        
        h = figure('DefaultAxesFontSize',16);
        h.Position(3) = round(h.Position(3)*1.2);
        bar(1:length(params.trigger_labels),[mean1; mean2]','grouped');
        hold on
        for k=1:length(params.trigger_labels)
            errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
            errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
        end
        
        p_values = zeros(n_triggers, 3);
        for i = 1:n_triggers
            [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
        end
        for i = 1:n_triggers
            sigstar({[i-0.11, i+0.11]}, p_values(i, 1));
        end
        
        hold off
        title('MNE: Group level peak power')
        ylabel('M60 pow')
        xlabel('Phalange')
        legend({'squidmag','opm'},'Location','eastoutside');
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_pow2_' peak_label cov hemi_labels{i_hemi} '.jpg']))
        
        %% Plot peak latencies
        data1 = lat_squidmag*1e3;
        data2 = lat_opm*1e3;
        data3 = lat_squidgrad*1e3;
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
        h.Position(3) = round(h.Position(3)*1.2);
        bar(1:length(params.trigger_labels),[mean1; mean2; mean3]','grouped');
        hold on
        for k=1:length(params.trigger_labels)
            errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
            errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
            errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
        end
        
        p_values = zeros(n_triggers, 3);
        for i = 1:n_triggers
            [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
            [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
            [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
        end
        for i = 1:n_triggers
            sigstar({[i-0.22, i]}, p_values(i, 1));
            sigstar({[i, i+0.22]}, p_values(i, 2));
            sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
        end
        
        hold off
        title('MNE: Group level peak latency')
        ylabel('M60 peak latency [ms]')
        xlabel('Phalange')
        legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
        xticklabels(params.trigger_labels)
        saveas(h, fullfile(base_save_path, 'figs', ['mne_latency_' peak_label cov hemi_labels{i_hemi} '.jpg']))
        
        data = {1e3*lat_squidmag, 1e3*lat_opm, 1e3*lat_squidgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak latency [ms]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Latency - SQMAG vs OPM vs SQGRAD'];
        save_path = fullfile(base_save_path, 'figs', ['mne_lat_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
        %% Distance between peak locations
        data = {dist_sqmag_opm, dist_sqgrad_opm, dist_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak distances [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak distances - SM-O vs SG-O vs SM-SG'];
        save_path = fullfile(base_save_path, 'figs', ['mne_dist_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);     
        close all    

        data = {distX_sqmag_opm, distX_sqgrad_opm, distX_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak distances [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak distances X - SM-O vs SG-O vs SM-SG'];
        save_path = fullfile(base_save_path, 'figs', ['mne_distX_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);     
        close all    

        data = {distY_sqmag_opm, distY_sqgrad_opm, distY_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak distances [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak distances Y - SM-O vs SG-O vs SM-SG'];
        save_path = fullfile(base_save_path, 'figs', ['mne_distY_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);     
        close all    

        data = {distZ_sqmag_opm, distZ_sqgrad_opm, distZ_sqmag_sqgrad};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Peak distances [mm]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak distances Z - SM-O vs SG-O vs SM-SG'];
        save_path = fullfile(base_save_path, 'figs', ['mne_distZ_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);     
        close all    
    
        %% Ratio of activation inside target region
        data = {targetregion_squidmag*1e2, targetregion_opm*1e2, targetregion_squidgrad*1e2};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Target region ratio [%]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak target region ratio'];
        save_path = fullfile(base_save_path, 'figs', ['mne_target_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);        
        close all

        %% Ratio of activation inside target region
        data = {overlap_opm_squidmag*1e2, overlap_opm_squidgrad*1e2, overlap_squidgrad_squidmag*1e2};
        triggerLabels = params.trigger_labels;
        yLabelStr = 'Overlap ratio [%]';
        titleStr = ['Group level ' params.peaks{1}.label ' MNE activation overlap - O-SM vs O-SG vs SG-SM'];
        save_path = fullfile(base_save_path, 'figs', ['mne_overlap_sqmag_opm_sqgrad' peak_label cov hemi_labels{i_hemi} '_box.jpg']);
        pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);        
        close all
    
        %% FAHM vs sub
%         for i_trigger = 1:n_triggers
%             h = figure('DefaultAxesFontSize',16);
%             plot(subs,fahm_squidmag(subs,i_trigger));
%             hold on
%             plot(subs,fahm_opm(subs,i_trigger));
%             plot(subs,fahm_squidgrad(subs,i_trigger));
%             hold off
%             title([params.trigger_labels{i_trigger} ' - mne FAHM over subjects'])
%             ylabel('FAHM [cm^2]')
%             xlabel('Subjects')
%             legend(['SQMAG '; 'OPM   '; 'SQGRAD'])
%             saveas(h, fullfile(base_save_path, 'figs', ['mne_fahm_vs_sub-' params.trigger_labels{i_trigger} '_' peak_label cov hemi_labels{i_hemi} '.jpg']))
%             close
%         end

        disp('done')
    end
end
end