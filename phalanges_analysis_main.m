%% Reset all
clear all
close all
restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    % Server:
    base_data_path = '/archive/21099_opm/';
    base_save_path = '/home/chrpfe/Documents/21099_opm/Phalanges';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/21099_opm/phalanges';
    on_server = true;
else
    % Laptop:
    base_data_path = '/Volumes/dataarchvie/21099_opm';
    base_save_path = '/Users/christophpfeiffer/data_local/Benchmarking_phalanges';
    base_matlab_path = '/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB';
    project_scripts_path = '/Users/christophpfeiffer/opm_phalanges';
    on_server = false;
end

%% Set up fieldtrip
addpath(fullfile(base_matlab_path,'fieldtrip/')) % Fieldtrip path
addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Overwrite
overwrite = [];
if on_server
    overwrite.preproc = true;
    overwrite.coreg = true;
    overwrite.mri = false;
    overwrite.dip = false;
    overwrite.empty_room = true;
    overwrite.mne = true;

    overwrite.sens_group = false;
    overwrite.dip_group = false;
    overwrite.mne_group = true;
else
    overwrite.preproc = false;
    overwrite.coreg = false;
    overwrite.mri = false;
    overwrite.dip = false;
    overwrite.empty_room = true;
    overwrite.mne = true;

    overwrite.sens_group = false;
    overwrite.dip_group = false;
    overwrite.mne_group = false;
end

%% Params
params = [];
params.pre = 0.03; %sec
params.post = 0.3; %sec
params.pad = 0.2; %sec

params.filter = [];
params.filter.hp_freq = 1;
params.filter.lp_freq = 70;
params.filter.bp_freq = [];
params.filter.notch = [50 60 100 120 150]; %[50 60 100 120 150];
params.do_hfc = false;
params.hfc_order = 1;
params.do_amm = true;
params.amm_in = 10;
params.amm_out = 2;
params.amm_thr = 1;

params.z_threshold = 20;
params.corr_threshold = 0.6; % correlation threshold for badchannel neighbors
params.opm_std_threshold = 5e-12;
params.eeg_std_threshold = 1e-4;
params.squidmag_std_threshold = 2.5e-12;
params.squidgrad_std_threshold = 5e-11;

params.n_comp = 40;
params.ica_cor = 0.8; % cutoff for EOG/ECG coherence
params.ica_coh = 0.95; % cutoff for EOG/ECG coherence

params.ds_freq = 500; % downsample frequency (timelock)

params.hpi_freq = 33;
params.hpi_gof = 0.9;

params.peaks = {};%cell(2,1);
params.peaks{1} = [];
params.peaks{1}.label = 'M60';
params.peaks{1}.peak_latency = [0.04 0.08];
params.peaks{2} = [];
params.peaks{2}.label = 'M100';
params.peaks{2}.peak_latency = [0.08 0.12];

params.trigger_code = [2 4 8 16 32];
params.phalange_labels = {'I3' 'I2' 'I1' 'T1' 'I2b'};

params.use_cov = 'empty_room'; % noise cov to use; default=prestim, alt: 'resting_state', 'all', 'empty_room'

%% Subjects + dates
subses = {'0005' '240208';
    '0905' '240229';
    '0916' '240320';
    '0953' '241104';
    '1096' '241022';
    '1153' '240321';
    '1167' '240425';
    '1186' '240925';
    '1190' '241023';
    '1191' '241024';
    '1193' '241029';
    '1194' '241029';
    '1195' '241030'};
mri_files = {'00000001.dcm' 
    '/mri/sub-15931_T1w.nii.gz'  
    '/nifti/anat/sub-15985_T1w.nii.gz'};

excl_subs = [];
if on_server
    subs_to_run = 1:size(subses,1);
else
    subs_to_run = 4; %1:size(subses,1)
end

%% Loop over subjects
for i_sub = setdiff(subs_to_run,excl_subs)
    params.sub = ['sub_' num2str(i_sub,'%02d')];

    %% Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    if ~exist(save_path, 'dir')
       mkdir(save_path)
    end
    if ~exist(fullfile(save_path,'figs'), 'dir')
       mkdir(fullfile(save_path,'figs'))
    end
    meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
    if i_sub == 9
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
    end
    opm_file = fullfile(raw_path, 'osmeg', 'PhalangesOPM_raw.fif');
    aux_file = fullfile(raw_path, 'meg', 'PhalangesEEG.fif');
    hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    %% OPM-MEG 
    if exist(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']),'file') && exist(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']),'file') && overwrite.preproc==false
        disp(['Not overwriting preproc for ' params.sub]);
    else
        ft_hastoolbox('mne', 1);

        if i_sub <=3 % Flip amplitudes in old recordings
            params.flip_sign  = true;
        else
            params.flip_sign  = false;
        end
        % Read data
        [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params); % Read data

        % OPM ICA
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        data_ica = ica_MEG(opm_cleaned, save_path, params, 1);
        clear opm_cleaned

        % OPM Average
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_ica

        % OPM-EEG ICA 
        cfg = [];
        cfg.elec = opmeeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_opmeeg_layout.mat']);
        opmeeg_layout = ft_prepare_layout(cfg);
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'opmeeg';
        data_ica = ica_MEG(opmeeg_cleaned, save_path, params, 1);
        close all
        clear opmeeg_cleaned

        params.modality = 'opmeeg';
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_ica

        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
        
        %% SQUID-MEG 
        ft_hastoolbox('mne', 1);

        % Read data
        [squid_cleaned, squideeg_cleaned] = read_cvMEG(meg_file, save_path, params); % Read data
        
        % SQUID-MAG ICA
        params.modality = 'squid';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'meg';
        data_ica = ica_MEG(squid_cleaned, save_path, params, 1);      

        % SQUID-MAG timelock
        params.modality = 'squidmag';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        timelock_MEG(data_ica, save_path, params);
        close all

        % SQUID-GRAD timelock
        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'megplanar';
        params.amp_scaler = 1e15/100;
        params.amp_label = 'B [fT/cm]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_ica

        % EEG ICA
        cfg = [];
        cfg.elec = squideeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_megeeg_layout.mat']);
        megeeg_layout = ft_prepare_layout(cfg);
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'squideeg';
        data_ica = ica_MEG(squideeg_cleaned, save_path, params, 1);
        close all
        clear squideeg_cleaned

        % EEG timelock
        params.modality = 'squideeg';
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        timelock_MEG(data_ica, save_path, params);
        close all
        clear data_íca

        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality    
    
        %% Save results in report
        create_bads_reports(base_save_path, i_sub, params);
    end
end

%% Sensor level group analysis
if overwrite.sens_group
    if ~exist(fullfile(base_save_path,'figs'), 'dir')
           mkdir(fullfile(base_save_path,'figs'))
    end
    subs = 1:13;
    sensor_results_goup(base_save_path,subs, params)
    close all
end

%% Exclude subs with missing co-reg
excl_subs = 1;

%% Prepare MRIs
for i_sub = setdiff(subs_to_run,excl_subs)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}],'mri');
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && overwrite.mri==false
        disp(['Not overwriting MRI for ' params.sub]);
    else
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        mri_file = fullfile(mri_path, 'orig','001.mgz');
        prepare_mri(mri_file,meg_file,save_path);
        close all
    end
end

%% HPI localization
for i_sub = setdiff(subs_to_run,excl_subs)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    hpi_path = fullfile(raw_path, 'osmeg'); %hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');

    if exist(fullfile(save_path, 'opm_trans.mat'),'file') && overwrite.coreg==false
        disp(['Not overwriting OPM transform for ' params.sub]);
    else
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        ft_hastoolbox('mne', 1);
        load(fullfile(save_path, [params.sub '_opm_ica_ds']));
        params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));
        fit_hpi(hpi_path, meg_file, save_path, params);
        close all
    end
end

%% Transform for OPM data
for i_sub = setdiff(subs_to_run,excl_subs)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);
    mri_path = fullfile(base_data_path,'MRI',['NatMEG_' subses{i_sub,1}]);
    if exist(fullfile(save_path, 'headmodels.mat'),'file') && exist(fullfile(save_path, 'opm_trans.mat'),'file') && or(overwrite.preproc,or(overwrite.coreg,overwrite.mri))
        clear headmodels meshes filename headshape 
        headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;
        meshes = load(fullfile(save_path,'meshes.mat')).meshes;    
        mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
        opm_trans = load(fullfile(save_path, 'opm_trans.mat')).opm_trans;
        
        clear opm_timelockedT opmeeg_timelcokedT 
        clear squideeg_timelocked squidmag_timelocked
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelocked.mat'])).timelocked;
        opmeeg_timelockedT = load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat'])).timelocked;
        squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat'])).timelocked;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked;

        % Transform opm & opmeeg data 
        meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', 'PhalangesMEG_proc-tsss+corr98.fif');
        end
        headshape = ft_read_headshape(meg_file);

        for i = 1:length(params.phalange_labels)
            opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.chanpos);
            opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.coilpos);
            opm_timelockedT{i}.grad.chanori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.chanori')';
            opm_timelockedT{i}.grad.coilori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.coilori')';
            opmeeg_timelockedT{i}.elec.chanpos = squideeg_timelocked{i}.elec.chanpos;
            opmeeg_timelockedT{i}.elec.elecpos = squideeg_timelocked{i}.elec.elecpos;
        end

        % Read and transform cortical restrained source model
        clear sourcemodel sourcemodel_inflated
        files = dir(fullfile(mri_path,'workbench'));
        if i_sub ==5
            files = dir(fullfile(save_path,'wb'));
        end
        for i = 1:length(files)
            if endsWith(files(i).name,'.L.midthickness.8k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            elseif endsWith(files(i).name,'.L.aparc.8k_fs_LR.label.gii')
                filename2 = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});

        aparc_L = ft_read_atlas({filename2,filename});
        aparc_R = ft_read_atlas({strrep(filename2,'.L.','.R.'),strrep(filename,'.L.','.R.')});
        tmp = ft_read_atlas(strrep(filename2, '.L.', '.R.'),'format','caret_label');
        n_labels = length(aparc_L.parcellationlabel);
        atlas = [];
        atlas.parcellationlabel = [aparc_L.parcellationlabel; aparc_R.parcellationlabel];
        atlas.parcellation = [aparc_L.parcellation; aparc_R.parcellation + n_labels];
        atlas.rgba = [aparc_L.rgba; aparc_R.rgba; [0 0 0 1]];
        n_labels = length(atlas.parcellationlabel);
        atlas.parcellation(isnan(atlas.parcellation))=n_labels+1;
        sourcemodel.brainstructure = atlas.parcellation;
        sourcemodel.brainstructurelabel = atlas.parcellationlabel;
        sourcemodel.brainstructurecolor = atlas.rgba;
        clear atlas aparc_L aparc_R

        T = mri_resliced.transform/mri_resliced.hdr.vox2ras;
        sourcemodel = ft_transform_geometry(T, sourcemodel);
        sourcemodel.inside = true(size(sourcemodel.pos,1),1);

        for i = 1:length(files)
            if endsWith(files(i).name,'.L.inflated.8k_fs_LR.surf.gii')
                filename = fullfile(mri_path,'workbench',files(i).name);
            end
        end
        sourcemodel_inflated = ft_read_headshape({filename, strrep(filename, '.L.', '.R.')});
        sourcemodel_inflated = ft_transform_geometry(T, sourcemodel_inflated);
        sourcemodel_inflated.inside = true(size(sourcemodel_inflated.pos,1),1);
        sourcemodel_inflated.brainstructure = sourcemodel.brainstructure;
        sourcemodel_inflated.brainstructurelabel = sourcemodel.brainstructurelabel;
        sourcemodel_inflated.brainstructurecolor = sourcemodel.brainstructurecolor;

        % Plot source and head models
        h=figure; 
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
        hold on; 
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
        ft_plot_sens(opmeeg_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
        hold off;
        title('OPM-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'figs', 'opm_layout2.fig'))
        saveas(h, fullfile(save_path, 'figs', 'opm_layout2.jpg'))
    
        h=figure; 
        ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
        hold on; 
        ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
        ft_plot_headmodel(headmodels.headmodel_meg, 'facealpha', 0.25, 'edgealpha', 0.25)
        ft_plot_sens(squidmag_timelocked{1}.grad,'unit','cm')
        ft_plot_sens(squideeg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
        hold off;
        title('SQUID-MEG')
        view([-140 10])
        savefig(h, fullfile(save_path, 'figs', 'meg_layout2.fig'))
        saveas(h, fullfile(save_path, 'figs', 'meg_layout2.jpg'))

        close all

        %% Save
        save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_opmeeg_timelockedT']), 'opmeeg_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_sourcemodel']), 'sourcemodel', '-v7.3');
        save(fullfile(save_path, [params.sub '_sourcemodel_inflated']), 'sourcemodel_inflated', '-v7.3');
        
        clear opm_timelockedT opmeeg_timelockedT sourcemodel_inflated sourcemodel
    else
        disp(['Required files not found. No transformed OPM/sourcemodel data was saved for ' params.sub])
    end
end

%% Dipole fits
for i_sub = setdiff(subs_to_run,excl_subs)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.sub);

    if exist(fullfile(save_path, 'dipoles.mat'),'file') && overwrite.dip==false
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
    else
        clear headmodels mri_resliced
        headmodels = load(fullfile(save_path, 'headmodels.mat')).headmodels;
        mri_resliced = load(fullfile(save_path, 'mri_resliced.mat')).mri_resliced;
        
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat'])).timelocked;
        
        for i_peak = 1:length(params.peaks)
            clear peak_opm peak_squidmag peak_squidgrad
            peak_opm = load(fullfile(save_path, [params.sub '_opm_' params.peaks{i_peak}.label])).peak; 
            peak_squidmag = load(fullfile(save_path, [params.sub '_squidmag_' params.peaks{i_peak}.label])).peak; 
            peak_squidgrad = load(fullfile(save_path, [params.sub '_squidgrad_' params.peaks{i_peak}.label])).peak; 
            fit_dipoles(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelockedT, headmodels, mri_resliced, peak_squidmag, peak_squidgrad, peak_opm, params);
            clear peak_opm peak_squidmag peak_squidgrad
        end
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
    end
end

%% Dipole group analysis
if overwrite.dip_group
    subs = setdiff(subs_to_run,excl_subs);
    dipole_results_goup(base_save_path,subs, params)
end

%% Empty room & resting state for noise covariances
for i_sub = setdiff(subs_to_run,excl_subs)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path,params.sub);
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    if i_sub <=3 % Flip amplitudes in old recordings
        params.flip_sign  = true;
    else
        params.flip_sign  = false;
    end
    if exist(fullfile(save_path, [params.sub '_resting_state_squid.mat']),'file') && overwrite.empty_room == false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    else
        clear data_ica
        data_ica = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT{1};
        opm_chs = data_ica.label(contains(data_ica.label,'bz'));
        opm_grad = data_ica.grad;
        clear data_ica
        data_ica = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked{1};
        squidmag_chs = data_ica.label(contains(data_ica.label,'MEG'));
        clear data_ica
        data_ica = load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat'])).timelocked{1};
        squidgrad_chs = data_ica.label(contains(data_ica.label,'MEG'));
        clear data_ica
        squid_chs = squidmag_chs;
        clear data_ica
        % Empty room
        opm_file = fullfile(raw_path, 'osmeg', 'EmptyRoomOPM_raw.fif');
        squid_file = fullfile(raw_path, 'meg', 'EmptyRoomMEG.fif');
        if exist(opm_file,'file') && exist(squid_file,'file')
            read_empty_rooms(opm_file, squid_file, opm_chs, squid_chs, opm_grad, save_path, params);
        end
    
        % RESO
        opm_file = fullfile(raw_path, 'osmeg', 'RSEOOPM_raw.fif');
        aux_file = fullfile(raw_path, 'meg', 'RSEOEEG.fif');
        read_osMEG_RS(opm_file, aux_file, opm_chs, save_path, params)
        squid_file = fullfile(raw_path, 'meg', 'RSEOMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
        read_cvMEG_RS(squid_file, squid_chs, save_path, params)
    end
end

%% Distributed source models
for i_sub = setdiff(subs_to_run,excl_subs)
    ft_hastoolbox('mne',1);
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path,params.sub);

    if exist(fullfile(save_path, 'opm_mne_peaks.mat'),'file') && overwrite.mne==false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    elseif i_sub == 1
        disp('SKIPPING SUBJECT - NO CO-REGISTRATION')
    else
        clear headmodels sourcemodel sourcemodel_inflated
        sourcemodel = load(fullfile(save_path, [params.sub '_sourcemodel'])).sourcemodel;
        sourcemodel_inflated = load(fullfile(save_path, [params.sub '_sourcemodel_inflated'])).sourcemodel_inflated;
        headmodels = load(fullfile(save_path,'headmodels.mat')).headmodels;

        %% Load data and cov
        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squidmag_timelocked = load(fullfile(save_path, [params.sub '_squidmag_timelocked.mat'])).timelocked;
        squidgrad_timelocked = load(fullfile(save_path, [params.sub '_squidgrad_timelocked.mat'])).timelocked;

        for i = 1:length(opm_timelockedT)
            opm_timelockedT{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_opm.mat'])).opm_RS_cov;
            squidmag_timelocked{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_squid.mat'])).squidmag_RS_cov;
            squidgrad_timelocked{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_squid.mat'])).squidgrad_RS_cov;   
            if exist(fullfile(save_path, [params.sub '_ER_squid.mat']),'file')
                opm_timelockedT{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_opm.mat'])).opm_ER_cov;
                squidmag_timelocked{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_squid.mat'])).squidmag_ER_cov;
                squidgrad_timelocked{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_squid.mat'])).squidgrad_ER_cov;
            end
        end
        
        %% MNE fit
        params.inv_method = 'mne';
        fit_mne(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelockedT, headmodels, sourcemodel, sourcemodel_inflated, params);

        %% ELORETA fit
        params.inv_method = 'eloreta';
        fit_eloreta(save_path, squidmag_timelocked, squidgrad_timelocked, opm_timelockedT, headmodels, sourcemodel, sourcemodel_inflated, params);

        clear squimdag_timelocked squidgrad_timelocked opm_timelockedT headmodels sourcemodel sourcemodel_inflated
    end
end
close all

%% MNE group analysis
if overwrite.mne_group
    subs = setdiff(subs_to_run,excl_subs);
    mne_results_goup(base_save_path, subs, params);
end

%%
close all
clear all
exit