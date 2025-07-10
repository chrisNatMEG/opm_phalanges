function mne_results_goup(base_save_path, subs, params)

cov = [];

for i_peak = 1:length(params.peaks)

    peak_label = params.peaks{i_peak}.label;

    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    n_subs = max(subs);
    n_ph = length(params.phalange_labels);
    dist_sqmag_opm = nan(n_subs,n_ph);
    dist_sqgrad_opm = nan(n_subs,n_ph);
    dist_sqmag_sqgrad = nan(n_subs,n_ph);
    fahm_opm = nan(n_subs,n_ph);
    fahm_squidmag = nan(n_subs,n_ph);
    fahm_squidgrad = nan(n_subs,n_ph);
    lat_opm = nan(n_subs,n_ph);
    lat_squidmag = nan(n_subs,n_ph);
    lat_squidgrad = nan(n_subs,n_ph);
    pow_opm = nan(n_subs,n_ph);
    pow_squidmag = nan(n_subs,n_ph);
    pow_squidgrad = nan(n_subs,n_ph);
    mom_opm = nan(n_subs,n_ph);
    mom_squidmag = nan(n_subs,n_ph);
    mom_squidgrad = nan(n_subs,n_ph);
    for i_sub = subs
        params.sub = ['sub_' num2str(i_sub,'%02d')];
        ft_hastoolbox('mne', 1);
        save_path = fullfile(base_save_path,params.sub);
        mne_squidmag{i_sub} = load(fullfile(save_path, 'squidmag_mne_peaks.mat')).squidmag_peak; 
        mne_squidgrad{i_sub} = load(fullfile(save_path, 'squidgrad_mne_peaks.mat')).squidgrad_peak;
        mne_opm{i_sub} = load(fullfile(save_path, 'opm_mne_peaks.mat')).opm_peak;
      
        % Metrics: 
        % - distance between mnes for same phalange different systems
        % - over phalanges: average distance from mean location within distance
    
        for i_phalange = 1:n_ph
            pos_squidmag(i_phalange,:) = mne_squidmag{i_sub}{i_phalange,i_peak}.loc;
            pos_squidgrad(i_phalange,:) = mne_squidgrad{i_sub}{i_phalange,i_peak}.loc;
            pos_opm(i_phalange,:) = mne_opm{i_sub}{i_phalange,i_peak}.loc;
    
            dist_sqmag_opm(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_opm(i_phalange,:));
            dist_sqgrad_opm(i_sub,i_phalange) = 1e1*norm(pos_squidgrad(i_phalange,:)-pos_opm(i_phalange,:));
            dist_sqmag_sqgrad(i_sub,i_phalange) = 1e1*norm(pos_squidmag(i_phalange,:)-pos_squidgrad(i_phalange,:));
    
            pow_squidmag(i_sub,i_phalange) = mne_squidmag{i_sub}{i_phalange,i_peak}.pow;
            pow_squidgrad(i_sub,i_phalange) = mne_squidgrad{i_sub}{i_phalange,i_peak}.pow;
            pow_opm(i_sub,i_phalange) = mne_opm{i_sub}{i_phalange,i_peak}.pow;
    
            mom_squidmag(i_sub,i_phalange) = mne_squidmag{i_sub}{i_phalange,i_peak}.mom;
            mom_squidgrad(i_sub,i_phalange) = mne_squidgrad{i_sub}{i_phalange,i_peak}.mom;
            mom_opm(i_sub,i_phalange) = mne_opm{i_sub}{i_phalange,i_peak}.mom;
    
            lat_squidmag(i_sub,i_phalange) = mne_squidmag{i_sub}{i_phalange,i_peak}.latency;
            lat_squidgrad(i_sub,i_phalange) = mne_squidgrad{i_sub}{i_phalange,i_peak}.latency;
            lat_opm(i_sub,i_phalange) = mne_opm{i_sub}{i_phalange,i_peak}.latency;
    
            fahm_opm(i_sub,i_phalange) = mne_opm{i_sub}{i_phalange,i_peak}.fahm; % mean distance from center of phalanges
            fahm_squidmag(i_sub,i_phalange) = mne_squidmag{i_sub}{i_phalange,i_peak}.fahm; % mean distance from center of phalanges
            fahm_squidgrad(i_sub,i_phalange) = mne_squidgrad{i_sub}{i_phalange,i_peak}.fahm; % mean distance from center of phalanges
    
            overlap_opm_squidmag(i_sub,i_phalange) = mne_opm{i_sub}{i_phalange,i_peak}.overlap_squidmag/mne_opm{i_sub}{i_phalange,i_peak}.fahm;
            overlap_opm_squidgrad(i_sub,i_phalange) = mne_opm{i_sub}{i_phalange,i_peak}.overlap_squidgrad/mne_opm{i_sub}{i_phalange,i_peak}.fahm;
            overlap_squidmag_squidgrad(i_sub,i_phalange) = mne_squidmag{i_sub}{i_phalange,i_peak}.overlap_squidgrad/mne_squidmag{i_sub}{i_phalange,i_peak}.fahm;
            overlap_squidmag_opm(i_sub,i_phalange) = mne_squidmag{i_sub}{i_phalange,i_peak}.overlap_opm/mne_squidmag{i_sub}{i_phalange,i_peak}.fahm;
            overlap_squidgrad_squidmag(i_sub,i_phalange) = mne_squidgrad{i_sub}{i_phalange,i_peak}.overlap_squidmag/mne_squidgrad{i_sub}{i_phalange,i_peak}.fahm;
            overlap_squidgrad_opm(i_sub,i_phalange) = mne_squidgrad{i_sub}{i_phalange,i_peak}.overlap_opm/mne_squidgrad{i_sub}{i_phalange,i_peak}.fahm;
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
    title(['Dist SQ-MAG to OPM (mean = ' num2str(mean(mean(dist_sqmag_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
    ylabel('Distance [mm]')
    xlabel('Phalange')
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_squidmag_to_opm_dist_' peak_label cov '.jpg']))
    
    h = figure('DefaultAxesFontSize',16);
    bar(1:length(params.phalange_labels),mean(dist_sqgrad_opm,1,'omitnan'));
    hold on
    er = errorbar(1:5,mean(dist_sqgrad_opm,1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-min(dist_sqgrad_opm,[],1,'omitnan'), mean(dist_sqgrad_opm,1,'omitnan')-max(dist_sqgrad_opm,[],1,'omitnan'));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    er.LineWidth = 1;
    er.CapSize = 30;
    hold off
    title(['Dist SQ-GRAD to OPM (mean = ' num2str(mean(mean(dist_sqgrad_opm,'omitnan'),'omitnan'),'%.1f') 'mm)'])
    ylabel('Distance [mm]')
    xlabel('Phalange')
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_squidgrad_to_opm_dist_' peak_label cov '.jpg']))
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
    title(['Dist SQ-MAG to SQ-GRAD (mean = ' num2str(mean(mean(dist_sqmag_sqgrad,'omitnan'),'omitnan'),'%.1f') 'mm)'])
    ylabel('Distance [mm]')
    xlabel('Phalange')
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_squidmag_to_squidgrad_dist_' peak_label cov '.jpg']))
    close
    
    for i_ph = 1:5
        h = figure('DefaultAxesFontSize',16);
        plot(subs,dist_sqmag_opm(subs,i_ph));
        hold on
        plot(subs,dist_sqgrad_opm(subs,i_ph));
        plot(subs,dist_sqmag_sqgrad(subs,i_ph));
        hold off
        title([params.phalange_labels{i_ph} ' - mne distances over subjects'])
        ylabel('Distance [mm]')
        xlabel('Subjects')
        legend(['SQMAG-OPM   '; 'SQGRAD-OPM  '; 'SQMAG-SQGRAD'])
        saveas(h, fullfile(base_save_path, 'figs', ['mne_dist_vs_sub-' params.phalange_labels{i_ph} '_' peak_label cov '.jpg']))
        close
    end
    
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
    bar(1:length(params.phalange_labels),[mean1; mean2; mean3]','grouped');
    hold on
    for k=1:length(params.phalange_labels)
        errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
        errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
    end
    
    p_values = zeros(5, 3);
    for i = 1:5
        [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
        [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
        [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
    end
    for i = 1:5
        sigstar({[i-0.22, i]}, p_values(i, 1));
        sigstar({[i, i+0.22]}, p_values(i, 2));
        sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
    end
    
    hold off
    title('MNE: Group level FAHM')
    ylabel('M60 FAHM [cm^2]')
    xlabel('Phalange')
    legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_fahm_' peak_label cov '.jpg']))
    
    data = {fahm_squidmag, fahm_opm, fahm_squidgrad};
    triggerLabels = params.phalange_labels;
    yLabelStr = 'FAHM [cm^2]';
    titleStr = ['Group level ' params.peaks{1}.label ' MNE FAHM - SQMAG vs OPM vs SQGRAD'];
    save_path = fullfile(base_save_path, 'figs', 'mne_fahm_sqmag_opm_sqgrad_box.jpg');
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
    bar(1:length(params.phalange_labels),[mean1; mean2; mean3]','grouped');
    hold on
    for k=1:length(params.phalange_labels)
        errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
        errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
    end
    
    p_values = zeros(5, 3);
    for i = 1:5
        [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
        [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
        [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
    end
    for i = 1:5
        sigstar({[i-0.22, i]}, p_values(i, 1));
        sigstar({[i, i+0.22]}, p_values(i, 2));
        sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
    end
    
    hold off
    title('MNE: Group level peak power')
    ylabel('M60 pow')
    xlabel('Phalange')
    legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_pow_' peak_label cov '.jpg']))
    
    data = {pow_squidmag, pow_opm, pow_squidgrad};
    triggerLabels = params.phalange_labels;
    yLabelStr = 'Peak power';
    titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Power - SQMAG vs OPM vs SQGRAD'];
    save_path = fullfile(base_save_path, 'figs', 'mne_pow_sqmag_opm_sqgrad_box.jpg');
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
    bar(1:length(params.phalange_labels),[mean1; mean2; mean3]','grouped');
    hold on
    for k=1:length(params.phalange_labels)
        errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
        errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
    end
    
    p_values = zeros(5, 3);
    for i = 1:5
        [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
        [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
        [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
    end
    for i = 1:5
        sigstar({[i-0.22, i]}, p_values(i, 1));
        sigstar({[i, i+0.22]}, p_values(i, 2));
        sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
    end
    
    hold off
    title('MNE: Group level peak moment')
    ylabel('M60 mom')
    xlabel('Phalange')
    legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_mom_' peak_label cov '.jpg']))
    
    data = {1e9*mom_squidmag, 1e9*mom_opm, 1e9*mom_squidgrad};
    triggerLabels = params.phalange_labels;
    yLabelStr = 'Peak moment [nAm]';
    titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Moment - SQMAG vs OPM vs SQGRAD'];
    save_path = fullfile(base_save_path, 'figs', 'mne_mom_sqmag_opm_sqgrad_box.jpg');
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);
    
    data = {1e9*mom_squidmag, 1e9*mom_opm};
    triggerLabels = params.phalange_labels;
    yLabelStr = 'Peak moment [nAm]';
    titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Moment - SQMAG vs OPM'];
    save_path = fullfile(base_save_path, 'figs', 'mne_mom_sqmag_opm_box.jpg');
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
    bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
    hold on
    for k=1:length(params.phalange_labels)
        errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    end
    
    p_values = zeros(5, 3);
    for i = 1:5
        [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
    end
    for i = 1:5
        sigstar({[i-0.11, i+0.11]}, p_values(i, 1));
    end
    
    hold off
    title('MNE: Group level peak power')
    ylabel('M60 pow')
    xlabel('Phalange')
    legend({'squidmag','opm'},'Location','eastoutside');
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_pow2_' peak_label cov '.jpg']))
    
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
    bar(1:length(params.phalange_labels),[mean1; mean2; mean3]','grouped');
    hold on
    for k=1:length(params.phalange_labels)
        errorbar(k-0.22,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
        errorbar(k,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
        errorbar(k+0.22,mean3(k),err3(1,k),err3(2,k),'k','linestyle','none');
    end
    
    p_values = zeros(5, 3);
    for i = 1:5
        [~, p_values(i, 1)] = ttest(data1(:,i), data2(:,i));
        [~, p_values(i, 2)] = ttest(data2(:,i), data3(:,i));
        [~, p_values(i, 3)] = ttest(data1(:,i), data3(:,i));
    end
    for i = 1:5
        sigstar({[i-0.22, i]}, p_values(i, 1));
        sigstar({[i, i+0.22]}, p_values(i, 2));
        sigstar({[i-0.22, i+0.22]}, p_values(i, 3));
    end
    
    hold off
    title('MNE: Group level peak latency')
    ylabel('M60 peak latency [ms]')
    xlabel('Phalange')
    legend({'squidmag','opm','squidgrad'},'Location','eastoutside');
    xticklabels(params.phalange_labels)
    saveas(h, fullfile(base_save_path, 'figs', ['mne_latency_' peak_label cov '.jpg']))
    
    data = {1e3*lat_squidmag, 1e3*lat_opm, 1e3*lat_squidgrad};
    triggerLabels = params.phalange_labels;
    yLabelStr = 'Peak latency [ms]';
    titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak Latency - SQMAG vs OPM vs SQGRAD'];
    save_path = fullfile(base_save_path, 'figs', 'mne_lat_sqmag_opm_sqgrad_box.jpg');
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

    %%
    data = {dist_sqmag_opm, dist_sqgrad_opm, dist_sqmag_sqgrad};
    triggerLabels = params.phalange_labels;
    yLabelStr = 'Peak distances [mm]';
    titleStr = ['Group level ' params.peaks{1}.label ' MNE Peak distances - SM-O vs SG-O vs SM-SG'];
    save_path = fullfile(base_save_path, 'figs', 'mne_dist_sqmag_opm_sqgrad_box.jpg');
    pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path,1);

    %% Plot FAHM squidgrad vs opm
    % data1 = fahm_squidgrad;
    % data2 = fahm_opm;
    % mean1 = mean(data1,1,'omitnan');
    % mean2 = mean(data2,1,'omitnan');
    % min1 = min(data1,[],1,'omitnan');
    % min2 = min(data2,[],1,'omitnan');
    % max1 = max(data1,[],1,'omitnan');
    % max2 = max(data2,[],1,'omitnan');
    % err1 = [mean1-min1; max1-mean1];
    % err2 = [mean2-min2; max2-mean2];
    % 
    % h = figure('DefaultAxesFontSize',16);
    % bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
    % hold on
    % for k=1:length(params.phalange_labels)
    %     errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    %     errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    % end
    % p_values = zeros(1, 5);
    % for i = 1:5
    %     [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
    % end
    % sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
    % hold off
    % title('MNE: Group level M100 FAHM')
    % ylabel('M60 FAHM [cm^2]')
    % xlabel('Phalange')
    % legend({'squidgrad','opm'});
    % xticklabels(params.phalange_labels)
    % saveas(h, fullfile(base_save_path, 'figs', 'mne_fahm_squidgrad_opm.jpg'))
    
    %% Plot FAHM squidgrad vs opm
    % data1 = fahm_squidgrad;
    % data2 = fahm_squidmag;
    % mean1 = mean(data1,1,'omitnan');
    % mean2 = mean(data2,1,'omitnan');
    % min1 = min(data1,[],1,'omitnan');
    % min2 = min(data2,[],1,'omitnan');
    % max1 = max(data1,[],1,'omitnan');
    % max2 = max(data2,[],1,'omitnan');
    % err1 = [mean1-min1; max1-mean1];
    % err2 = [mean2-min2; max2-mean2];
    % 
    % h = figure('DefaultAxesFontSize',16);
    % bar(1:length(params.phalange_labels),[mean1; mean2]','grouped');
    % hold on
    % for k=1:length(params.phalange_labels)
    %     errorbar(k-0.15,mean1(k),err1(1,k),err1(2,k),'k','linestyle','none');
    %     errorbar(k+0.15,mean2(k),err2(1,k),err2(2,k),'k','linestyle','none');
    % end
    % p_values = zeros(1, 5);
    % for i = 1:5
    %     [~, p_values(i)] = ttest(data1(:, i), data2(:, i));
    % end
    % sigstar({[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]}, p_values);
    % hold off
    % title('MNE: Group level M100 FAHM')
    % ylabel('M60 FAHM [cm^2]')
    % xlabel('Phalange')
    % legend({'squidgrad','squidmag'});
    % xticklabels(params.phalange_labels)
    % saveas(h, fullfile(base_save_path, 'figs', 'mne_fahm_squidgrad_squidmag.jpg'))
    
    close all
    
    %% FAHM vs sub
    for i_ph = 1:5
        h = figure('DefaultAxesFontSize',16);
        plot(subs,fahm_squidmag(subs,i_ph));
        hold on
        plot(subs,fahm_opm(subs,i_ph));
        plot(subs,fahm_squidgrad(subs,i_ph));
        hold off
        title([params.phalange_labels{i_ph} ' - mne FAHM over subjects'])
        ylabel('FAHM [cm^2]')
        xlabel('Subjects')
        legend(['SQMAG '; 'OPM   '; 'SQGRAD'])
        saveas(h, fullfile(base_save_path, 'figs', ['mne_fahm_vs_sub-' params.phalange_labels{i_ph} '_' peak_label cov '.jpg']))
        close
    end
    disp('done')
end
end