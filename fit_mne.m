function fit_mne(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelocked, headmodels, sourcemodel, sourcemodel_inflated, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Prepare leadfields
headmodel = headmodels.headmodel_meg;
sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';
sourcemodel.unit = 'cm';

%% Load peaks
squidmag_peaks = cell(length(params.peaks),1);
squidgrad_peaks = cell(length(params.peaks),1);
opm_peaks = cell(length(params.peaks),1);
% for i_peak = 1:length(params.peaks)
%     squidmag_peaks{i_peak} = load(fullfile(save_path, [params.sub '_squidmag_' params.peaks{i_peak}.label])).M60;%peak; 
%     squidgrad_peaks{i_peak} = load(fullfile(save_path, [params.sub '_squidgrad_' params.peaks{i_peak}.label])).M60;%peak; 
%     opm_peaks{i_peak} = load(fullfile(save_path, [params.sub '_opm_' params.peaks{i_peak}.label])).M60;%peak; 
%     for i_ph = 1:5
%         squidmag_peaks{i_peak}{i_ph}.label = 'M60';
%         squidgrad_peaks{i_peak}{i_ph}.label = 'M60';
%         opm_peaks{i_peak}{i_ph}.label = 'M60';
%     end
% end

%% MNE invserse
squidmag_mne_peak = cell(length(params.trigger_code),length(params.peaks));
squidgrad_mne_peak = cell(length(params.trigger_code),length(params.peaks));
opm_mne_peak = cell(length(params.trigger_code),length(params.peaks));

for i_phalange = 1:length(params.trigger_code)
    params.i_phalange = i_phalange;
    cov = '';
    if isfield(params,'use_cov') && strcmp(params.use_cov,'all')
        squidmag_timelocked{i_phalange}.cov = squidmag_timelocked{i_phalange}.cov_all;
        squidgrad_timelocked{i_phalange}.cov = squidgrad_timelocked{i_phalange}.cov_all;
        opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_all;
        cov = '_covAll';
    elseif isfield(params,'use_cov') && strcmp(params.use_cov,'resting_state')
        squidmag_timelocked{i_phalange}.cov = squidmag_timelocked{i_phalange}.cov_RS;
        squidgrad_timelocked{i_phalange}.cov = squidgrad_timelocked{i_phalange}.cov_RS;
        opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_RS;
        cov = '_covRS';
    elseif isfield(params,'use_cov') && strcmp(params.use_cov,'empty_room')
        if ~isfield(squidmag_timelocked{i_phalange},'cov_ER')
            continue % skip if no empty room covariance available
        end
        squidmag_timelocked{i_phalange}.cov = squidmag_timelocked{i_phalange}.cov_ER;
        squidgrad_timelocked{i_phalange}.cov = squidgrad_timelocked{i_phalange}.cov_ER;
        opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_ER;
        cov = '_covER';
    end

    %% MEG-MAG
    cfg = [];
    cfg.grad             = squidmag_timelocked{i_phalange}.grad; % sensor positions
    cfg.senstype         = 'meg';            % sensor type
    cfg.sourcemodel      = sourcemodel;           % source points
    cfg.headmodel        = headmodel;          % volume conduction model
    cfg.normalize        = 'yes';
    leadfield = ft_prepare_leadfield(cfg,squidmag_timelocked{i_phalange});

    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    tmp = ft_sourceanalysis(cfg, squidmag_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;

    params.modality = 'squidmag';

    for i_peak = 1:length(params.peaks)
        squidmag_mne_peak{i_phalange,i_peak} = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak},params,save_path);
    
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = squidmag_mne_peak{i_phalange,i_peak}.latency;
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, tmp)
        lighting gouraud
        material dull
        title(['SQUID-MAG (FAHM=' num2str(squidmag_mne_peak{i_phalange,i_peak}.fahm,3) 'cm^2; t=' num2str(round(squidmag_mne_peak{i_phalange,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_' squidmag_mne_peak{i_phalange,i_peak}.label '_mne_ph' params.phalange_labels{i_phalange} cov '.jpg']))
        close all
    end

    clear tmp leadfield

    %% MEG-GRAD
    cfg = [];
    cfg.grad             = squidgrad_timelocked{i_phalange}.grad;              % sensor positions
    cfg.senstype         = 'meg';            % sensor type
    cfg.sourcemodel      = sourcemodel;           % source points
    cfg.headmodel        = headmodel;          % volume conduction model
    cfg.normalize        = 'yes';
    leadfield = ft_prepare_leadfield(cfg,squidgrad_timelocked{i_phalange});

    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.sourcemodel         = leadfield;
    cfg.keepfilter          = 'yes';
    tmp = ft_sourceanalysis(cfg, squidgrad_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;

    params.modality = 'squidgrad';

    for i_peak = 1:length(params.peaks)
        squidgrad_mne_peak{i_phalange,i_peak} = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak},params,save_path);
    
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = squidgrad_mne_peak{i_phalange,i_peak}.latency;
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, tmp)
        lighting gouraud
        material dull
        title(['SQUID-GRAD (FAHM=' num2str(squidgrad_mne_peak{i_phalange,i_peak}.fahm,3) 'cm^2; t=' num2str(round(squidgrad_mne_peak{i_phalange,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_' squidgrad_mne_peak{i_phalange,i_peak}.label '_mne_ph' params.phalange_labels{i_phalange} cov '.jpg']))
        close all
    end

    clear tmp leadfield

    %% OPM
    cfg = [];
    cfg.grad             = opm_timelocked{i_phalange}.grad;              % sensor positions
    cfg.senstype         = 'meg';            % sensor type
    cfg.sourcemodel      = sourcemodel;           % source points
    cfg.headmodel        = headmodel;          % volume conduction model
    cfg.normalize        = 'yes';
    leadfield = ft_prepare_leadfield(cfg,opm_timelocked{i_phalange});

    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;

    params.modality = 'opm';
    
    for i_peak = 1:length(params.peaks)
        opm_mne_peak{i_phalange,i_peak} = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak},params,save_path);
    
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = opm_mne_peak{i_phalange,i_peak}.latency;
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, tmp)
        lighting gouraud
        material dull
        title(['OPM (FAHM=' num2str(opm_mne_peak{i_phalange,i_peak}.fahm,3) 'cm^2; t=' num2str(round(opm_mne_peak{i_phalange,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_opm_' opm_mne_peak{i_phalange,i_peak}.label '_mne_ph' params.phalange_labels{i_phalange} cov '.jpg']))
        close all
    end

    clear tmp leadfield

    %% Overlaps
    for i_peak = 1:length(params.peaks)
        i_vertices = opm_mne_peak{i_phalange,i_peak}.halfmax_distribution & squidmag_mne_peak{i_phalange,i_peak}.halfmax_distribution;
        [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
        triangles = sourcemodel.tri(triangles,:);
        opm_mne_peak{i_phalange,i_peak}.overlap_squidmag = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
        squidmag_mne_peak{i_phalange,i_peak}.overlap_opm = opm_mne_peak{i_phalange,i_peak}.overlap_squidmag;
    
        i_vertices = opm_mne_peak{i_phalange,i_peak}.halfmax_distribution & squidgrad_mne_peak{i_phalange,i_peak}.halfmax_distribution;
        [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
        triangles = sourcemodel.tri(triangles,:);
        opm_mne_peak{i_phalange,i_peak}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
        squidgrad_mne_peak{i_phalange,i_peak}.overlap_opm = opm_mne_peak{i_phalange,i_peak}.overlap_squidgrad;
    
        i_vertices = squidmag_mne_peak{i_phalange,i_peak}.halfmax_distribution & squidgrad_mne_peak{i_phalange,i_peak}.halfmax_distribution;
        [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
        triangles = sourcemodel.tri(triangles,:);
        squidmag_mne_peak{i_phalange,i_peak}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
        squidgrad_mne_peak{i_phalange,i_peak}.overlap_squidmag = squidmag_mne_peak{i_phalange,i_peak}.overlap_squidgrad;
    end

end
save(fullfile(save_path, ['squidmag_mne_peaks' cov]), 'squidmag_mne_peak'); 
save(fullfile(save_path, ['squidgrad_mne_peaks' cov]), 'squidgrad_mne_peak'); 
save(fullfile(save_path, ['opm_mne_peaks' cov]), 'opm_mne_peak'); 

end