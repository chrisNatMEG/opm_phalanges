function fit_mne(save_path, squid_timelocked, opm_timelocked, headmodel, sourcemodel, sourcemodel_inflated, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(params,'plot_inflated')
    params.plot_inflated = false;
end
n_triggers = length(opm_timelocked);

headmodel = ft_convert_units(headmodel,'cm');
sourcemodel = ft_convert_units(sourcemodel,'cm');
if params.source_fixedori
    sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';
end

%% Prepare leadfields
cfg = [];
cfg.grad             = squid_timelocked{1}.grad; % sensor positions
cfg.channel          = 'megmag';
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
%cfg.normalize        = 'yes';
leadfield_squidmag = ft_prepare_leadfield(cfg,squid_timelocked{1});

cfg = [];
cfg.grad             = squid_timelocked{1}.grad; % sensor positions
cfg.channel          = 'meggrad';
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
%cfg.normalize        = 'yes';
leadfield_squidgrad = ft_prepare_leadfield(cfg,squid_timelocked{1});

cfg = [];
cfg.grad             = opm_timelocked{1}.grad; % sensor positions
cfg.senstype         = 'meg';            % sensor type
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
%cfg.normalize        = 'yes';
leadfield_opm = ft_prepare_leadfield(cfg,opm_timelocked{1});

%% Loop over triggers
squidmag_peak = cell(length(n_triggers),length(params.peaks));
squidgrad_peak = cell(length(n_triggers),length(params.peaks));
opm_peak = cell(length(n_triggers),length(params.peaks));

for i_trigger = 1:n_triggers
    params.i_trigger = i_trigger;

    %% Set covariance matrix (based on params.noise_cov selection)
    cov = '';
    if isfield(params,'noise_cov')
        if strcmp(params.noise_cov,'all')
            cov = '_covAll';
            squid_timelocked{i_trigger}.cov = squid_timelocked{i_trigger}.cov_all;
            opm_timelocked{i_trigger}.cov = opm_timelocked{i_trigger}.cov_all;
        elseif strcmp(params.noise_cov,'resting_state') && ~isfield(squid_timelocked{i_trigger},'cov_RS')
            cov = '_covRS';
            if isfield(squid_timelocked{i_trigger},'cov_RS') && size(opm_timelocked{i_trigger}.cov_RS,1) == size(opm_timelocked{i_trigger}.cov,1)
                squid_timelocked{i_trigger}.cov = squid_timelocked{i_trigger}.cov_RS;
                opm_timelocked{i_trigger}.cov = opm_timelocked{i_trigger}.cov_RS;
            else
                warning('Resting state covariance not existing or incorrect size.');
                break
            end
        elseif strcmp(params.noise_cov,'empty_room') && ~isfield(squid_timelocked{i_trigger},'cov_ER')
            cov = '_covER';
            if isfield(squid_timelocked{i_trigger},'cov_ER') && size(opm_timelocked{i_trigger}.cov_ER) == size(opm_timelocked{i_trigger}.cov)
                squid_timelocked{i_trigger}.cov = squid_timelocked{i_trigger}.cov_ER;
                opm_timelocked{i_trigger}.cov = opm_timelocked{i_trigger}.cov_ER;
            else
                warning('Empty room covariance not existing or incorrect size.');
                break
            end
        end
    end

    %% MEG-MAG
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_squidmag;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = 'megmag';
    tmp = ft_sourceanalysis(cfg, squid_timelocked{i_trigger});
    tmp.tri = sourcemodel.tri;

    params.modality = 'squidmag';

    if iscell(tmp.avg.mom) && size(tmp.avg.mom{1},1)~=3
        if size(tmp.avg.mom,1) == 1
            tmp.avg.mom = cell2mat(tmp.avg.mom');
        else
            tmp.avg.mom = cell2mat(tmp.avg.mom);
        end
    end
    
    h = figure;
    plot(tmp.time*1e3,tmp.avg.pow)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.trigger_labels{params.i_trigger}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.inv_method '_sourcepow_trig-' params.trigger_labels{params.i_trigger} '.jpg']))
    close all

    if params.plot_inflated
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
    end

    for i_peak = 1:length(params.peaks)
        peak = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak}.peak_latency,params);
        peak.label = params.peaks{i_peak}.label;

        h = plot_source_distribution(tmp, peak, params); 
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' peak.label '_mne_trig-' params.trigger_labels{i_trigger} cov '.jpg']))
        close all 

        squidmag_peak{i_trigger,i_peak} = peak;
        clear peak
    end
    clear tmp

    %% MEG-GRAD
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.sourcemodel         = leadfield_squidgrad;
    cfg.keepfilter          = 'yes';
    cfg.channel             = 'meggrad';
    tmp = ft_sourceanalysis(cfg, squid_timelocked{i_trigger});
    tmp.tri = sourcemodel.tri;

    params.modality = 'squidgrad';

    if iscell(tmp.avg.mom) && size(tmp.avg.mom{1},1)~=3
        if size(tmp.avg.mom,1) == 1
            tmp.avg.mom = cell2mat(tmp.avg.mom');
        else
            tmp.avg.mom = cell2mat(tmp.avg.mom);
        end
    end
    h = figure;
    plot(tmp.time*1e3,tmp.avg.pow)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.trigger_labels{params.i_trigger}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.inv_method '_sourcepow_trig-' params.trigger_labels{params.i_trigger} '.jpg']))
    close all

    if params.plot_inflated
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
    end
    
    for i_peak = 1:length(params.peaks)
        peak = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak}.peak_latency,params);
        peak.label = params.peaks{i_peak}.label;

        h = plot_source_distribution(tmp, peak, params); 
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' peak.label '_mne_trig-' params.trigger_labels{i_trigger} cov '.jpg']))
        close all  

        squidgrad_peak{i_trigger,i_peak} = peak;
        clear peak
    end
    clear tmp

    %% OPM
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_opm;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = '*bz';
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_trigger});
    tmp.tri = sourcemodel.tri;

    params.modality = 'opm';

    if iscell(tmp.avg.mom) && size(tmp.avg.mom{1},1)~=3
        if size(tmp.avg.mom,1) == 1
            tmp.avg.mom = cell2mat(tmp.avg.mom');
        else
            tmp.avg.mom = cell2mat(tmp.avg.mom);
        end
    end
    h = figure;
    plot(tmp.time*1e3,tmp.avg.pow)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.trigger_labels{params.i_trigger}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.inv_method '_sourcepow_trig-' params.trigger_labels{params.i_trigger} '.jpg']))
    close all
    
    if params.plot_inflated
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
    end
    
    for i_peak = 1:length(params.peaks)
        peak = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak}.peak_latency,params);
        peak.label = params.peaks{i_peak}.label;

        h = plot_source_distribution(tmp, peak, params); 
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' peak.label '_mne_trig-' params.trigger_labels{i_trigger} cov '.jpg']))
        close all  

        opm_peak{i_trigger,i_peak} = peak;
        clear peak
    end
    clear tmp


    %% Overlaps
    if size(opm_peak{i_trigger,i_peak}.loc,1) == 2
        for i_peak = 1:length(params.peaks)
            for i_hemi = 1:2
                i_vertices = intersect(opm_peak{i_trigger,i_peak}.halfmax_distribution{i_hemi},squidmag_peak{i_trigger,i_peak}.halfmax_distribution{i_hemi});
                [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
                triangles = sourcemodel.tri(triangles,:);
                opm_peak{i_trigger,i_peak}.overlap_squidmag(i_hemi) = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
                squidmag_peak{i_trigger,i_peak}.overlap_opm(i_hemi) = opm_peak{i_trigger,i_peak}.overlap_squidmag(i_hemi);
            
                i_vertices = intersect(opm_peak{i_trigger,i_peak}.halfmax_distribution{i_hemi}, squidgrad_peak{i_trigger,i_peak}.halfmax_distribution{i_hemi});
                [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
                triangles = sourcemodel.tri(triangles,:);
                opm_peak{i_trigger,i_peak}.overlap_squidgrad(i_hemi) = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
                squidgrad_peak{i_trigger,i_peak}.overlap_opm(i_hemi) = opm_peak{i_trigger,i_peak}.overlap_squidgrad(i_hemi);
            
                i_vertices = intersect(squidmag_peak{i_trigger,i_peak}.halfmax_distribution{i_hemi}, squidgrad_peak{i_trigger,i_peak}.halfmax_distribution{i_hemi});
                [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
                triangles = sourcemodel.tri(triangles,:);
                squidmag_peak{i_trigger,i_peak}.overlap_squidgrad(i_hemi) = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
                squidgrad_peak{i_trigger,i_peak}.overlap_squidmag(i_hemi) = squidmag_peak{i_trigger,i_peak}.overlap_squidgrad(i_hemi);
            end
        end
    else
        for i_peak = 1:length(params.peaks)
            i_vertices = intersect(opm_peak{i_trigger,i_peak}.halfmax_distribution,squidmag_peak{i_trigger,i_peak}.halfmax_distribution);
            [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
            triangles = sourcemodel.tri(triangles,:);
            opm_peak{i_trigger,i_peak}.overlap_squidmag = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
            squidmag_peak{i_trigger,i_peak}.overlap_opm = opm_peak{i_trigger,i_peak}.overlap_squidmag;
        
            i_vertices = intersect(opm_peak{i_trigger,i_peak}.halfmax_distribution, squidgrad_peak{i_trigger,i_peak}.halfmax_distribution);
            [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
            triangles = sourcemodel.tri(triangles,:);
            opm_peak{i_trigger,i_peak}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
            squidgrad_peak{i_trigger,i_peak}.overlap_opm = opm_peak{i_trigger,i_peak}.overlap_squidgrad;
        
            i_vertices = intersect(squidmag_peak{i_trigger,i_peak}.halfmax_distribution, squidgrad_peak{i_trigger,i_peak}.halfmax_distribution);
            [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
            triangles = sourcemodel.tri(triangles,:);
            squidmag_peak{i_trigger,i_peak}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
            squidgrad_peak{i_trigger,i_peak}.overlap_squidmag = squidmag_peak{i_trigger,i_peak}.overlap_squidgrad;
        end
    end

end

if ~isempty(squidmag_peak{1,1})
    save(fullfile(save_path, ['squidmag_mne_peaks' cov]), 'squidmag_peak'); 
    save(fullfile(save_path, ['squidgrad_mne_peaks' cov]), 'squidgrad_peak'); 
    save(fullfile(save_path, ['opm_mne_peaks' cov]), 'opm_peak'); 
end

end