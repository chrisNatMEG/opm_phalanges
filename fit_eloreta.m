function fit_eloreta(save_path, squid_timelocked, opm_timelocked, headmodel, sourcemodel, sourcemodel_inflated, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isfield(params,'plot_inflated')
    params.plot_inflated = false;
end

%% Prepare leadfields
headmodel = ft_convert_units(headmodel,'cm');
sourcemodel = ft_convert_units(sourcemodel,'cm');
if params.source_fixedori
    sourcemodel.mom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex')';
end

%% invserse
squidmag_peak = cell(length(params.trigger_code),length(params.peaks));
squidgrad_peak = cell(length(params.trigger_code),length(params.peaks));
opm_peak = cell(length(params.trigger_code),length(params.peaks));

%% Leadfields
cfg = [];
cfg.grad             = squid_timelocked{1}.grad; % sensor positions
cfg.senstype         = 'meg';            % sensor type
cfg.channel          = 'megmag';
cfg.sourcemodel      = sourcemodel;           % source points
cfg.headmodel        = headmodel;          % volume conduction model
%cfg.normalize        = 'yes';
leadfield_squidmag = ft_prepare_leadfield(cfg,squid_timelocked{1});

cfg = [];
cfg.grad             = squid_timelocked{1}.grad; % sensor positions
cfg.senstype         = 'meg';            % sensor type
cfg.channel          = 'meggrad';
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

for i_phalange = 1:length(params.trigger_code)
    params.i_phalange = i_phalange;
    cov = '';
    if isfield(params,'use_cov') && strcmp(params.use_cov,'all')
        squid_timelocked{i_phalange}.cov = squid_timelocked{i_phalange}.cov_all;
        opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_all;
        cov = '_covAll';
    elseif isfield(params,'use_cov') && strcmp(params.use_cov,'resting_state')
        if size(opm_timelocked{i_phalange}.cov_RS) == size(opm_timelocked{i_phalange}.cov)
            squid_timelocked{i_phalange}.cov = squid_timelocked{i_phalange}.cov_RS;
            opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_RS;
            cov = '_covRS';
        else
            continue
        end
    elseif isfield(params,'use_cov') && strcmp(params.use_cov,'empty_room')
        if ~isfield(squid_timelocked{i_phalange},'cov_ER')
            continue % skip if no empty room covariance available
        end
        squid_timelocked{i_phalange}.cov = squid_timelocked{i_phalange}.cov_ER;
        opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_ER;
        cov = '_covER';
        if size(opm_timelocked{i_phalange}.cov_ER) == size(opm_timelocked{i_phalange}.cov)
            squid_timelocked{i_phalange}.cov = squid_timelocked{i_phalange}.cov_ER;
            opm_timelocked{i_phalange}.cov = opm_timelocked{i_phalange}.cov_ER;
            cov = '_covER';
        else 
            continue
        end
    end

    %% MEG-MAG
    cfg = [];
    cfg.method              = 'eloreta';
    cfg.eloreta.prewhiten       = 'yes';
    cfg.eloreta.lambda          = 3;
    cfg.eloreta.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_squidmag;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = 'megmag';
    tmp = ft_sourceanalysis(cfg, squid_timelocked{i_phalange});
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
    plot(tmp.time*1e3,tmp.avg.mom.^2)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.phalange_labels{params.i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.inv_method '_sourcepow_ph-' params.phalange_labels{params.i_phalange} '.jpg']))
    close all

    for i_peak = 1:length(params.peaks)
        squidmag_peak{i_phalange,i_peak} = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak},params,save_path);
    end
    if params.plot_inflated
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
    end
    for i_peak = 1:length(params.peaks)
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = squidmag_peak{i_phalange,i_peak}.latency;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, tmp)
        lighting gouraud
        material dull
        title(['SQUID-MAG (FAHM=' num2str(squidmag_peak{i_phalange,i_peak}.fahm,3) 'cm^2; t=' num2str(round(squidmag_peak{i_phalange,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_' squidmag_peak{i_phalange,i_peak}.label '_eloreta_ph' params.phalange_labels{i_phalange} cov '.jpg']))
        close all
    end

    clear tmp

    %% MEG-GRAD
    cfg = [];
    cfg.method              = 'eloreta';
    cfg.eloreta.prewhiten       = 'yes';
    cfg.eloreta.lambda          = 3;
    cfg.eloreta.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.sourcemodel         = leadfield_squidgrad;
    cfg.keepfilter          = 'yes';
    cfg.channel             = 'meggrad';
    tmp = ft_sourceanalysis(cfg, squid_timelocked{i_phalange});
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
    plot(tmp.time*1e3,tmp.avg.mom.^2)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.phalange_labels{params.i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.inv_method '_sourcepow_ph-' params.phalange_labels{params.i_phalange} '.jpg']))
    close all

    for i_peak = 1:length(params.peaks)
        squidgrad_peak{i_phalange,i_peak} = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak},params,save_path);
    end
    if params.plot_inflated
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
    end
    for i_peak = 1:length(params.peaks)
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = squidgrad_peak{i_phalange,i_peak}.latency;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, tmp)
        lighting gouraud
        material dull
        title(['SQUID-GRAD (FAHM=' num2str(squidgrad_peak{i_phalange,i_peak}.fahm,3) 'cm^2; t=' num2str(round(squidgrad_peak{i_phalange,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_' squidgrad_peak{i_phalange,i_peak}.label '_eloreta_ph' params.phalange_labels{i_phalange} cov '.jpg']))
        close all
    end

    clear tmp

    %% OPM
    cfg = [];
    cfg.method              = 'eloreta';
    cfg.eloreta.prewhiten       = 'yes';
    cfg.eloreta.lambda          = 3;
    cfg.eloreta.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodel;    % supply the headmodel
    cfg.sourcemodel         = leadfield_opm;
    cfg.senstype            = 'meg';            % sensor type
    cfg.keepfilter          = 'yes';
    cfg.channel             = '*bz';
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
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
    plot(tmp.time*1e3,tmp.avg.mom.^2)
    xlabel('t [msec]')
    ylabel('Field power')
    xlim([-params.pre params.post]*1e3);
    title([params.modality ' - ' params.phalange_labels{params.i_phalange}])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' params.inv_method '_sourcepow_ph-' params.phalange_labels{params.i_phalange} '.jpg']))
    close all
    
    for i_peak = 1:length(params.peaks)
        opm_peak{i_phalange,i_peak} = FullAreaHalfMax(tmp,sourcemodel,params.peaks{i_peak},params,save_path);
    end
    if params.plot_inflated
        tmp.pos = sourcemodel_inflated.pos;
        tmp.tri = sourcemodel_inflated.tri;
    end
    for i_peak = 1:length(params.peaks)
        cfg = [];
        cfg.method          = 'surface';
        cfg.funparameter    = 'pow';
        cfg.funcolormap     = 'jet';    
        cfg.colorbar        = 'no';
        cfg.latency         = opm_peak{i_phalange,i_peak}.latency;
        h = figure;
        h.Position(3) = round(h.Position(3)*1.2);
        ft_sourceplot(cfg, tmp)
        lighting gouraud
        material dull
        title(['OPM (FAHM=' num2str(opm_peak{i_phalange,i_peak}.fahm,3) 'cm^2; t=' num2str(round(opm_peak{i_phalange,i_peak}.latency*1e3)) 'ms)'])
        saveas(h, fullfile(save_path,'figs', [params.sub '_opm_' opm_peak{i_phalange,i_peak}.label '_eloreta_ph' params.phalange_labels{i_phalange} cov '.jpg']))
        close all
    end

    clear tmp

    %% Overlaps
    for i_peak = 1:length(params.peaks)
        i_vertices = opm_peak{i_phalange,i_peak}.halfmax_distribution & squidmag_peak{i_phalange,i_peak}.halfmax_distribution;
        [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
        triangles = sourcemodel.tri(triangles,:);
        opm_peak{i_phalange,i_peak}.overlap_squidmag = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
        squidmag_peak{i_phalange,i_peak}.overlap_opm = opm_peak{i_phalange,i_peak}.overlap_squidmag;
    
        i_vertices = opm_peak{i_phalange,i_peak}.halfmax_distribution & squidgrad_peak{i_phalange,i_peak}.halfmax_distribution;
        [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
        triangles = sourcemodel.tri(triangles,:);
        opm_peak{i_phalange,i_peak}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
        squidgrad_peak{i_phalange,i_peak}.overlap_opm = opm_peak{i_phalange,i_peak}.overlap_squidgrad;
    
        i_vertices = squidmag_peak{i_phalange,i_peak}.halfmax_distribution & squidgrad_peak{i_phalange,i_peak}.halfmax_distribution;
        [triangles,~] = find(ismember(sourcemodel.tri,i_vertices)); 
        triangles = sourcemodel.tri(triangles,:);
        squidmag_peak{i_phalange,i_peak}.overlap_squidgrad = sum(calculateTriangleAreas(sourcemodel.pos, triangles))/3;
        squidgrad_peak{i_phalange,i_peak}.overlap_squidmag = squidmag_peak{i_phalange,i_peak}.overlap_squidgrad;
    end

end
save(fullfile(save_path, ['squidmag_eloreta_peaks' cov]), 'squidmag_peak'); 
save(fullfile(save_path, ['squidgrad_eloreta_peaks' cov]), 'squidgrad_peak'); 
save(fullfile(save_path, ['opm_eloreta_peaks' cov]), 'opm_peak'); 

end