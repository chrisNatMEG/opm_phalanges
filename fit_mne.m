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
        if strcmp(params.noise_cov,'empty_room') && ~isfield(squid_timelocked{i_trigger},'cov_ER')
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
    %cfg.latency            = [0 params.post];
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
    plot(tmp.time*1e3,std(tmp.avg.pow,0,1))
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

        h = plot_source_distribution(tmp, peak, params,0); 
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' peak.label '_mne_trig-' params.trigger_labels{i_trigger} cov '.jpg']))
        close all 

        squidmag_peak{i_trigger,i_peak} = peak;
        clear peak
    end
    squidmag_mne{i_trigger} = tmp;
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
    %cfg.latency            = [0 params.post];
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
    plot(tmp.time*1e3,std(tmp.avg.pow,0,1))
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

        h = plot_source_distribution(tmp, peak, params, 0); 
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' peak.label '_mne_trig-' params.trigger_labels{i_trigger} cov '.jpg']))
        close all  

        squidgrad_peak{i_trigger,i_peak} = peak;
        clear peak
    end
    squidgrad_mne{i_trigger} = tmp;
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
    %cfg.latency            = [0 params.post];
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
    plot(tmp.time*1e3,std(tmp.avg.pow,0,1))
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

        h = plot_source_distribution(tmp, peak, params, 0); 
        saveas(h, fullfile(save_path,'figs', [params.sub '_' params.modality '_' peak.label '_mne_trig-' params.trigger_labels{i_trigger} cov '.jpg']))
        close all  

        opm_peak{i_trigger,i_peak} = peak;
        clear peak
    end
    opm_mne{i_trigger} = tmp;
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

colors = [[0 0.4470 0.7410]; % blue
    [0.8500 0.3250 0.0980]; % red
    [0.9290 0.6940 0.1250]; % yellow
    [0.4940 0.1840 0.5560]; % purple
    [0.4660 0.6740 0.1880]; % green
    [0.6350 0.0780 0.1840]]; % light blue

%% SQMAG
params.modality = 'squidmag';
peak_label = params.peaks{1}.label;
tmp = squidmag_mne{1};
tmp.avg.pow = tmp.avg.pow/n_triggers;
peak = squidmag_peak{1,1};
peak.latency  = peak.latency/n_triggers;
peak.fahm = peak.fahm/n_triggers;
for i_trigger = 2:n_triggers
    tmp.avg.pow = tmp.avg.pow + squidmag_mne{i_trigger}.avg.pow/n_triggers;
    peak.latency = peak.latency + squidmag_peak{i_trigger,1}.latency/n_triggers;
    peak.fahm = peak.fahm + squidmag_peak{i_trigger,1}.fahm/n_triggers;
end
h = plot_source_distribution(tmp, peak, params, 0);
subplot(1,2,1)
hold on
for i_trigger = 1:n_triggers
    for i_pk = 1:length(squidmag_peak{i_trigger,1}.iloc)
        ft_plot_dipole(tmp.pos(squidmag_peak{i_trigger,1}.iloc(i_pk),:),[0 0 0],'color',colors(i_trigger,:))
    end
end
hold off
subplot(1,2,2)
hold on
for i_trigger = 1:n_triggers
    for i_pk = 1:length(squidmag_peak{i_trigger,1}.iloc)
        ft_plot_dipole(tmp.pos(squidmag_peak{i_trigger,1}.iloc(i_pk),:),[0 0 0],'color',colors(i_trigger,:))
    end
end
hold off
saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_' peak_label '_mne' cov '.jpg']))
close all 

%% SQGRAD
params.modality = 'squidgrad';
tmp = squidgrad_mne{1};
tmp.avg.pow = tmp.avg.pow/n_triggers;
peak = squidgrad_peak{1,1};
peak.latency  = peak.latency/n_triggers;
peak.fahm = peak.fahm/n_triggers;
for i_trigger = 2:n_triggers
    tmp.avg.pow = tmp.avg.pow + squidgrad_mne{i_trigger}.avg.pow/n_triggers;
    peak.latency = peak.latency + squidgrad_peak{i_trigger,1}.latency/n_triggers;
    peak.fahm = peak.fahm + squidgrad_peak{i_trigger,1}.fahm/n_triggers;
end
h = plot_source_distribution(tmp, peak, params, 0);
subplot(1,2,1)
hold on
for i_trigger = 1:n_triggers
    for i_pk = 1:length(squidgrad_peak{i_trigger,1}.iloc)
        ft_plot_dipole(tmp.pos(squidgrad_peak{i_trigger,1}.iloc(i_pk),:),[0 0 0],'color',colors(i_trigger,:))
    end
end
subplot(1,2,2)
hold on
for i_trigger = 1:n_triggers
    for i_pk = 1:length(squidgrad_peak{i_trigger,1}.iloc)
        ft_plot_dipole(tmp.pos(squidgrad_peak{i_trigger,1}.iloc(i_pk),:),[0 0 0],'color',colors(i_trigger,:))
    end
end
hold off
saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_' peak_label '_mne' cov '.jpg']))
close all

%% OPM
params.modality = 'opm';
tmp = opm_mne{1};
tmp.avg.pow = tmp.avg.pow/n_triggers;
peak = opm_peak{1,1};
peak.latency  = peak.latency/n_triggers;
peak.fahm = peak.fahm/n_triggers;
for i_trigger = 2:n_triggers
    tmp.avg.pow = tmp.avg.pow + opm_mne{i_trigger}.avg.pow/n_triggers;
    peak.latency = peak.latency + opm_peak{i_trigger,1}.latency/n_triggers;
    peak.fahm = peak.fahm + opm_peak{i_trigger,1}.fahm/n_triggers;
end
h = plot_source_distribution(tmp, peak, params, 0);
subplot(1,2,1)
hold on
for i_trigger = 1:n_triggers
    for i_pk = 1:length(opm_peak{i_trigger,1}.iloc)
        ft_plot_dipole(tmp.pos(opm_peak{i_trigger,1}.iloc(i_pk),:),[0 0 0],'color',colors(i_trigger,:))
    end
end
hold off
subplot(1,2,2)
hold on
for i_trigger = 1:n_triggers
    for i_pk = 1:length(opm_peak{i_trigger,1}.iloc)
        ft_plot_dipole(tmp.pos(opm_peak{i_trigger,1}.iloc(i_pk),:),[0 0 0],'color',colors(i_trigger,:))
    end
end
hold off
saveas(h, fullfile(save_path,'figs', [params.sub '_opm_' peak_label '_mne' cov '.jpg']))
close all

%% Save
save(fullfile(save_path, ['mne_distributions' cov]), 'squidmag_mne','squidgrad_mne','opm_mne'); 

if ~isempty(squidmag_peak{1,1})
    peaks = squidmag_peak;
    save(fullfile(save_path, ['squidmag_mne_peaks' cov]), 'peaks'); 
    peaks = squidgrad_peak;
    save(fullfile(save_path, ['squidgrad_mne_peaks' cov]), 'peaks'); 
    peaks = opm_peak;
    save(fullfile(save_path, ['opm_mne_peaks' cov]), 'peaks'); 
end

if isfield(params,'do_sourcemovie') && params.do_sourcemovie
    % OPM
    for i_trigger = 1:n_triggers
        v = VideoWriter(fullfile(save_path,'figs', [params.sub '_opm_' params.trigger_labels{i_trigger} '.avi']),"Uncompressed AVI");
        v.FrameRate = 10;
        
        v.open
        tmp = opm_mne{i_trigger};
        maxpow = max(max(tmp.avg.pow));
        [~,i_start] = min(abs(tmp.time-0));
        for i = i_start:length(tmp.time)
            
            viewangles = [90 0 -90 0];
            cfg = [];
            cfg.method          = 'surface';
            cfg.funparameter    = 'pow';
            cfg.funcolormap     = 'jet';    
            cfg.colorbar        = 'no';
            cfg.latency         = tmp.time(i);
            h = figure;
            subplot(1,2,1); % right hemisphere
            cfg.figure = h;
            ft_sourceplot(cfg, tmp)
            clim([0 maxpow])
            material dull
            view(viewangles(1),viewangles(2))
            camlight();
            lighting gouraud
            title([' (t=' num2str(1e3*tmp.time(i),3) 'ms)'])
            subplot(1,2,2); % left hemisphere
            cfg.figure = h;
            ft_sourceplot(cfg, tmp)
            clim([0 maxpow])
            material dull
            view(viewangles(3),viewangles(4))
            camlight()
            lighting gouraud
            title([' (t=' num2str(1e3*tmp.time(i),3) 'ms)'])
            set(h,'Position',[10 10 970 450]);
            writeVideo(v,getframe(h));
            close all
        end
        v.close
    end
    
    % SQ-GRAD
    for i_trigger = 1:n_triggers
        v = VideoWriter(fullfile(save_path,'figs', [params.sub '_sqgrad_' params.trigger_labels{i_trigger} '.avi']),"Uncompressed AVI");
        v.FrameRate = 10;
    
        v.open
        tmp = squidgrad_mne{i_trigger};
        maxpow = max(max(tmp.avg.pow));
        [~,i_start] = min(abs(tmp.time-0));
        for i = i_start:length(tmp.time)
            
            viewangles = [90 0 -90 0];
            cfg = [];
            cfg.method          = 'surface';
            cfg.funparameter    = 'pow';
            cfg.funcolormap     = 'jet';    
            cfg.colorbar        = 'no';
            cfg.latency         = tmp.time(i);
            h = figure;
            subplot(1,2,1); % right hemisphere
            cfg.figure = h;
            ft_sourceplot(cfg, tmp)
            clim([0 maxpow])
            material dull
            view(viewangles(1),viewangles(2))
            camlight();
            lighting gouraud
            title([' (t=' num2str(1e3*tmp.time(i),3) 'ms)'])
            subplot(1,2,2); % left hemisphere
            cfg.figure = h;
            ft_sourceplot(cfg, tmp)
            clim([0 maxpow])
            material dull
            view(viewangles(3),viewangles(4))
            camlight()
            lighting gouraud
            title([' (t=' num2str(1e3*tmp.time(i),3) 'ms)'])
            set(h,'Position',[10 10 970 450]);
            writeVideo(v,getframe(h));
            close all
        end
        v.close
    end

    % SQ-MAG
    for i_trigger = 1:n_triggers
        v = VideoWriter(fullfile(save_path,'figs', [params.sub '_sqmag_' params.trigger_labels{i_trigger} '.avi']),"Uncompressed AVI");
        v.FrameRate = 10;
    
        v.open
        tmp = squidmag_mne{i_trigger};
        maxpow = max(max(tmp.avg.pow));
        [~,i_start] = min(abs(tmp.time-0));
        for i = i_start:length(tmp.time)
            
            viewangles = [90 0 -90 0];
            cfg = [];
            cfg.method          = 'surface';
            cfg.funparameter    = 'pow';
            cfg.funcolormap     = 'jet';    
            cfg.colorbar        = 'no';
            cfg.latency         = tmp.time(i);
            h = figure;
            subplot(1,2,1); % right hemisphere
            cfg.figure = h;
            ft_sourceplot(cfg, tmp)
            clim([0 maxpow])
            material dull
            view(viewangles(1),viewangles(2))
            camlight();
            lighting gouraud
            title([' (t=' num2str(1e3*tmp.time(i),3) 'ms)'])
            subplot(1,2,2); % left hemisphere
            cfg.figure = h;
            ft_sourceplot(cfg, tmp)
            clim([0 maxpow])
            material dull
            view(viewangles(3),viewangles(4))
            camlight()
            lighting gouraud
            title([' (t=' num2str(1e3*tmp.time(i),3) 'ms)'])
            set(h,'Position',[10 10 970 450]);
            writeVideo(v,getframe(h));
            close all
        end
        v.close
    end
end

end