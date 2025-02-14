function fit_mne(save_path,squidmag_timelocked,squidgrad_timelocked,squideeg_timelocked,opm_timelocked,opmeeg_timelocked,headmodels,sourcemodel,latency,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Prepare leadfields
headmodel = headmodels.headmodel_meg;
headmodel.order = 20;
cfg = [];
cfg.grad             = squidmag_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'squidmag';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodel;          % volume conduction model
leadfield_squidmag = ft_prepare_leadfield(cfg,squidmag_timelocked{1});

cfg = [];
cfg.grad             = squidgrad_timelocked{1}.grad;              % sensor positions
cfg.channel          = 'squidgrad';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodel;          % volume conduction model
leadfield_squidgrad = ft_prepare_leadfield(cfg,squidgrad_timelocked{1});

cfg = [];
cfg.grad             = opm_timelocked{1}.grad;              % sensor positions
cfg.channel          = '*bz';                  % the used channels
cfg.senstype         = 'meg';            % sensor type
cfg.grid.pos         = sourcemodel.pos;           % source points
cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
cfg.headmodel        = headmodel;          % volume conduction model
leadfield_opm = ft_prepare_leadfield(cfg,opm_timelocked{1});

% if ~isempty(headmodels.headmodel_eeg)
%     cfg = [];
%     cfg.elec             = squideeg_timelocked{1}.elec;              % sensor positions
%     cfg.channel          = 'eeg';                  % the used channels
%     cfg.senstype         = 'eeg';            % sensor type
%     cfg.grid.pos         = sourcemodel.pos;           % source points
%     cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
%     cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
%     leadfield_squideeg = ft_prepare_leadfield(cfg,squideeg_timelocked{1});
%     
%     cfg = [];
%     cfg.elec             = opmeeg_timelocked{1}.elec;              % sensor positions
%     cfg.channel          = 'eeg';                  % the used channels
%     cfg.senstype         = 'eeg';            % sensor type
%     cfg.grid.pos         = sourcemodel.pos;           % source points
%     cfg.grid.inside      = sourcemodel.inside; % all source points are inside of the brain
%     cfg.headmodel        = headmodels.headmodel_eeg;          % volume conduction model
%     leadfield_opmeeg = ft_prepare_leadfield(cfg,opmeeg_timelocked{1});
% end

%% MNE invserse
% MEG-MAG
squidmag_mne = [];
squidmag_mne.avg = cell(5,1);
squidmag_mne_M100 = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodels.headmodel_meg;    % supply the headmodel
    cfg.sourcemodel         = leadfield_squidmag;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = 'megmag';         % which channels to use
    tmp = ft_sourceanalysis(cfg, squidmag_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    cfg = [];
    cfg.projectmom = 'yes';
    tmp = ft_sourcedescriptives(cfg,tmp);
    squidmag_mne.avg{i_phalange} = [];
    squidmag_mne.avg{i_phalange}.pow = tmp.avg.pow;
    squidmag_mne.avg{i_phalange}.mom = tmp.avg.mom;
    squidmag_mne_M100{i_phalange} = [];
    [squidmag_mne_M100{i_phalange}.fahm, squidmag_mne_M100{i_phalange}.peakloc] = FullAreaHalfMax(tmp,sourcemodel,latency{i_phalange}.squidmag);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = latency{i_phalange}.squidmag;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['SQUID-MAG (FAHM=' num2str(squidmag_mne_M100{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidmag_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
squidmag_mne.time = tmp.time;
squidmag_mne.cfg = tmp.cfg;
squidmag_mne.method = tmp.method;
squidmag_mne.pos = tmp.pos;
squidmag_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'squidmag_mne'), 'squidmag_mne'); 
save(fullfile(save_path, 'squidmag_mne_M100'), 'squidmag_mne_M100'); 
clear tmp squidmag_mne leadfield_squidmag

% MEG-GRAD
squidgrad_mne = [];
squidgrad_mne.avg = cell(5,1);
squidgrad_mne_M100 = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodels.headmodel_meg;    % supply the headmodel
    cfg.senstype            = 'meg';
    cfg.channel             = 'megplanar';            % which channels to use
    cfg.sourcemodel         = leadfield_squidgrad;
    tmp = ft_sourceanalysis(cfg, squidgrad_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    cfg = [];
    cfg.projectmom = 'yes';
    tmp = ft_sourcedescriptives(cfg,tmp);
    squidgrad_mne.avg{i_phalange} = [];
    squidgrad_mne.avg{i_phalange}.pow = tmp.avg.pow;
    squidgrad_mne.avg{i_phalange}.mom = tmp.avg.mom;
    squidgrad_mne_M100{i_phalange} = [];
    [squidgrad_mne_M100{i_phalange}.fahm, squidgrad_mne_M100{i_phalange}.peakloc] = FullAreaHalfMax(tmp,sourcemodel,latency{i_phalange}.squidgrad);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = latency{i_phalange}.squidgrad;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['SQUID-GRAD (FAHM=' num2str(squidgrad_mne_M100{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_squidgrad_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
squidgrad_mne.time = tmp.time;
squidgrad_mne.cfg = tmp.cfg;
squidgrad_mne.method = tmp.method;
squidgrad_mne.pos = tmp.pos;
squidgrad_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'squidgrad_mne'), 'squidgrad_mne'); 
save(fullfile(save_path, 'squidgrad_mne_M100'), 'squidgrad_mne_M100'); 
clear tmp squidgrad_mne leadfield_squidgrad

% % MEG-EEG
% squideeg_mne = [];
% squideeg_mne.avg = cell(5,1);
% squideeg_mne_M100 = cell(5,1);
% for i_phalange = 1:5
%     squideeg_mne.avg{i_phalange} = [];
%     squideeg_mne_M100{i_phalange} = [];
%     if ~isempty(headmodels.headmodel_eeg)
%         try
%             cfg = [];
%             cfg.method              = 'mne';
%             cfg.mne.prewhiten       = 'yes';
%             cfg.mne.lambda          = 3;
%             cfg.mne.scalesourcecov  = 'yes';
%             cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
%             cfg.senstype            = 'eeg';            % sensor type
%             cfg.channel             = 'eeg';         % which channels to use
%             %cfg.sourcemodel.leaddield = leadfield_squideeg.leadfield;
%             cfg.sourcemodel         = leadfield_squideeg;
%             tmp = ft_sourceanalysis(cfg, squideeg_timelocked{i_phalange});
%             tmp.tri = sourcemodel.tri;
%             squideeg_mne.avg{i_phalange}.pow = tmp.avg.pow;
%             squideeg_mne.avg{i_phalange}.mom = tmp.avg.mom;
%             squideeg_mne.time = tmp.time;
%             squideeg_mne.cfg = tmp.cfg;
%             squideeg_mne.method = tmp.method;
%             squideeg_mne.pos = tmp.pos;
%             squideeg_mne.tri = sourcemodel.tri;
%             squideeg_mne_M100{i_phalange} = [];
%             [squideeg_mne_M100{i_phalange}.fahm, squideeg_mne_M100{i_phalange}.peakloc] = FullAreaHalfMax(tmp,sourcemodel,latency{i_phalange}.squideeg);
%     
%             cfg = [];
%             cfg.method          = 'surface';
%             cfg.funparameter    = 'pow';
%             cfg.funcolormap     = 'jet';    
%             cfg.colorbar        = 'no';
%             cfg.latency         = latency{i_phalange}.squideeg;
%             h = figure;
%             ft_sourceplot(cfg, tmp)
%             title(['SQUID-EEG (FAHM=' num2str(squideeg_mne_M100{i_phalange}.fahm,3) ')'])
%             saveas(h, fullfile(save_path,'figs', [params.sub '_squideeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
%             close all
%         catch
%             disp(['ERROR in SQUID-EEG ' params.sub ' - ph' num2str(i_phalange)])
%         end
%     end
% end
% 
% save(fullfile(save_path, 'squideeg_mne'), 'squideeg_mne'); 
% save(fullfile(save_path, 'squideeg_mne_M100'), 'squideeg_mne_M100'); 
% clear tmp squideeg_mne leadfield_squideeg

% OPM
opm_mne = [];
opm_mne.avg = cell(5,1);
opm_mne_M100 = cell(5,1);
for i_phalange = 1:5
    cfg = [];
    cfg.method              = 'mne';
    cfg.mne.prewhiten       = 'yes';
    cfg.mne.lambda          = 3;
    cfg.mne.scalesourcecov  = 'yes';
    cfg.headmodel           = headmodels.headmodel_meg;    % supply the headmodel
    cfg.sourcemodel         = leadfield_opm;
    cfg.senstype            = 'meg';            % sensor type
    cfg.channel             = '*bz';         % which channels to use
    tmp = ft_sourceanalysis(cfg, opm_timelocked{i_phalange});
    tmp.tri = sourcemodel.tri;
    cfg = [];
    cfg.projectmom = 'yes';
    tmp = ft_sourcedescriptives(cfg,tmp);
    opm_mne.avg{i_phalange} = [];
    opm_mne.avg{i_phalange}.pow = tmp.avg.pow;
    opm_mne.avg{i_phalange}.mom = tmp.avg.mom;
    opm_mne_M100{i_phalange} = [];
    [opm_mne_M100{i_phalange}.fahm, opm_mne_M100{i_phalange}.peakloc] = FullAreaHalfMax(tmp,sourcemodel,latency{i_phalange}.opm);

    cfg = [];
    cfg.method          = 'surface';
    cfg.funparameter    = 'pow';
    cfg.funcolormap     = 'jet';    
    cfg.colorbar        = 'no';
    cfg.latency         = latency{i_phalange}.opm;
    h = figure;
    ft_sourceplot(cfg, tmp)
    title(['OPM (FAHM=' num2str(opm_mne_M100{i_phalange}.fahm,3) ')'])
    saveas(h, fullfile(save_path,'figs', [params.sub '_opm_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
    close all
end
opm_mne.time = tmp.time;
opm_mne.cfg = tmp.cfg;
opm_mne.method = tmp.method;
opm_mne.pos = tmp.pos;
opm_mne.tri = sourcemodel.tri;

save(fullfile(save_path, 'opm_mne'), 'opm_mne'); 
save(fullfile(save_path, 'opm_mne_M100'), 'opm_mne_M100'); 
clear tmp opm_mne leadfield_opm

% % OPM-EEG
% opmeeg_mne = [];
% opmeeg_mne.avg = cell(5,1);
% opmeeg_mne_M100 = cell(5,1);
% for i_phalange = 1:5
%     opmeeg_mne.avg{i_phalange} = [];
%     opmeeg_mne_M100{i_phalange} = [];
%     if ~isempty(headmodels.headmodel_eeg)
%         try
%             cfg = [];
%             cfg.method              = 'mne';
%             cfg.mne.prewhiten       = 'yes';
%             cfg.mne.lambda          = 3;
%             cfg.mne.scalesourcecov  = 'yes';
%             cfg.headmodel           = headmodels.headmodel_eeg;    % supply the headmodel
%             cfg.senstype            = 'eeg';            % sensor type
%             cfg.channel             = 'eeg';         % which channels to use
%             cfg.sourcemodel         = leadfield_opmeeg;
%             tmp = ft_sourceanalysis(cfg, opmeeg_timelocked{i_phalange});
%             tmp.tri = sourcemodel.tri;
%             opmeeg_mne.avg{i_phalange}.pow = tmp.avg.pow;
%             opmeeg_mne.avg{i_phalange}.mom = tmp.avg.mom;
%             opmeeg_mne.time = tmp.time;
%             opmeeg_mne.cfg = tmp.cfg;
%             opmeeg_mne.method = tmp.method;
%             opmeeg_mne.pos = tmp.pos;
%             opmeeg_mne.tri = sourcemodel.tri;
%             opmeeg_mne_M100{i_phalange} = [];
%             [opmeeg_mne_M100{i_phalange}.fahm, opmeeg_mne_M100{i_phalange}.peakloc] = FullAreaHalfMax(tmp,sourcemodel,latency{i_phalange}.opmeeg);
%             
%             cfg = [];
%             cfg.method          = 'surface';
%             cfg.funparameter    = 'pow';
%             cfg.funcolormap     = 'jet';    
%             cfg.colorbar        = 'no';
%             cfg.latency         = latency{i_phalange}.opmeeg;
%             h = figure;
%             ft_sourceplot(cfg, tmp)
%             title(['OPM-EEG (FAHM=' num2str(opmeeg_mne_M100{i_phalange}.fahm,3) ')'])
%             saveas(h, fullfile(save_path,'figs', [params.sub '_opmeeg_mne_ph' params.phalange_labels{i_phalange} '.jpg']))
%             close all
%         catch
%             disp(['ERROR in OPM-EEG ' params.sub ' - ph' num2str(i_phalange)])
%         end
%     end
% end
% 
% save(fullfile(save_path, 'opmeeg_mne'), 'opmeeg_mne'); 
% save(fullfile(save_path, 'opmeeg_mne_M100'), 'opmeeg_mne_M100'); 
% clear tmp opmeeg_mne leadfield_opmeeg

end