%% Reset all
clear all
close all
restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    % Server:
    base_data_path = '/archive/21099_opm/';
    base_save_path = '/home/chrpfe/Documents/21099_opm/';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/21099_opm/phalanges';
    on_server = true;
else
    % Laptop:
    base_data_path = '/Volumes/dataarchvie/21099_opm';
    base_save_path = '/Users/christophpfeiffer/data_local/Benchmarking/';
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
    overwrite.dip = true;
    overwrite.empty_room = true;
    overwrite.mne = true;

    overwrite.sens_group = true;
    overwrite.dip_group = true;
    overwrite.mne_group = true;
else
    overwrite.preproc = false;
    overwrite.coreg = false;
    overwrite.mri = false;
    overwrite.dip = true;
    overwrite.empty_room = false;
    overwrite.mne = false;

    overwrite.sens_group = false;
    overwrite.dip_group = false;
    overwrite.mne_group = false;
end

%% Params
params = [];
params.paradigm = 'Phalanges';

% Trialse
params.pre = 0.03; %sec
params.post = 0.3; %sec
params.pad = 0.2; %sec
params.delay = 0.041; % Stimulus delay in seconds (e.g., 0.01 for eartubes or 0.041 for membranes).

% Filter
params.filter = [];
params.filter.hp_freq = 0.1;
params.filter.lp_freq = 70;
params.filter.bp_freq = [];
params.filter.notch = [50 60 100 120 150]; %[50 60 100 120 150];

% Spatiotemporal filter (OPM-MEG only)
params.do_hfc = true;
params.hfc_order = 1;
params.do_amm = false;
params.amm_in = 12;
params.amm_out = 2;
params.amm_thr = 1;

% Bad channel and trial detection thresholds
params.outlier_zscore = 3; % Outliers: how many stddevs above mean
params.outlier_ratio = 0.5; % Outliers: ratio of frequencies over threshold above which the channel is considered an outlier
params.corr_threshold = 0.6; % Correlation threshold for badchannei detection based on neighbor correlation
params.z_threshold = 20; % Zmax threshold for badchannel and trial detection based on jumps
% params.opm_std_threshold = 5e-12; % Standard deviation threshold for OPM badtrial detection
% params.squidmag_std_threshold = 2.5e-12; % Standard deviation for SQUID-MAG badtrial detection
% params.squidgrad_std_threshold = 2000e-13; % Standard deviation for SQUID-GRAD badtrial detection
% params.eeg_std_threshold = 10e-4; % Standard deviation for EEG badtrial detection
params.opm_range_threshold = 20e-12; % Range for OPM badtrial detection
params.squidmag_range_threshold = 10e-12; % Range for SQUID-MAG badtrial detection
params.squidgrad_range_threshold = 4000e-13; % Range for SQUID-GRAD badtrial detection


% ICA ECG&EOG artifact removal 
params.n_comp = 40;
params.ica_cor = 0.8; % cutoff for EOG/ECG coherence
params.ica_coh = 0.95; % cutoff for EOG/ECG coherence

% Timelocking
params.ds_freq = 500; % downsample frequency (timelock)
params.trigger_codes = {2 4 8 16 32};
params.trigger_labels = {'I3' 'I2' 'I1' 'T1' 'I2b'};
params.peaks = {};
params.peaks{1} = [];
params.peaks{1}.label = 'M60';
params.peaks{1}.peak_latency = [0.04 0.08];

% HPI coregistration
params.hpi_freq = 33;
params.hpi_gof = 0.9;

% Source reconstruction - dipoles
params.numdipoles = 1;

% Source reconstruction - distributed
params.source_fixedori = true; 
params.covs = {' '};%{'empty_room', ' '}; % noise cov to use; default=prestim, alt: 'resting_state', 'all', 'empty_room' , prestim = ' '
params.mne_view = 'sides';
params.plot_inflated = true;
params.target_region = 'postcentral';

save(fullfile(base_save_path, params.paradigm), 'params', '-v7.3');          

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
    '1195' '241030';
    '1209' '250219';
    '1215' '250415'};

bads =  {[];
    [];
    [];
    [];
    [];
    {'R403_bz', 'R408_bz', 'R409_bz'};
    {'R403_bz', 'R408_bz', 'R409_bz'};
    {'L209_bz', 'R209_bz', 'R403_bz', 'R408_bz'};
    {'L209_bz', 'R408_bz'};
    [];
    [];
    [];
    [];
    [];
    []};

if on_server
    subs_to_run = 1:size(subses,1);
else
    subs_to_run = 2; %1:size(subses,1)
end
excl_subs = [];
excl_subs_src = excl_subs;

%% Loop over subjects
for i_sub = setdiff(subs_to_run,excl_subs)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    params.manual_bads = bads{i_sub}';
    if i_sub <=3 % Flip amplitudes in old recordings
        params.flip_sign  = true;
    else
        params.flip_sign  = false;
    end

    %% Paths
    raw_path = fullfile(base_data_path,'MEG',['NatMEG_' subses{i_sub,1}], subses{i_sub,2});
    save_path = fullfile(base_save_path,params.paradigm,params.sub);
    save_path_mri = fullfile(base_save_path,'MRI',params.sub);
    hpi_path = fullfile(raw_path, 'osmeg'); %hpi_file = fullfile(raw_path, 'osmeg', 'HPIpre_raw.fif');
    
    % Create folders if they do not exist yet
    if ~exist(fullfile(base_save_path,params.paradigm), 'dir')
        mkdir(fullfile(base_save_path,params.paradigm))
    end
    if ~exist(save_path, 'dir')
       mkdir(save_path)
    end
    if ~exist(fullfile(save_path,'figs'), 'dir')
       mkdir(fullfile(save_path,'figs'))
    end

    meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif']);
    if ~exist(meg_file,'file')
        meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_proc-tsss+corr98.fif']);
    end
    opm_file = fullfile(raw_path, 'osmeg', [params.paradigm 'OPM_raw.fif']);
    aux_file = fullfile(raw_path, 'meg', [params.paradigm 'EEG.fif']);
    
    %% OPM-MEG 
    if exist(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat']),'file') && exist(fullfile(save_path, [params.sub '_squideeg_timelocked.mat']),'file') && overwrite.preproc==false
        disp(['Not overwriting preproc for ' params.sub]);
    else
        ft_hastoolbox('mne', 1);

        % Read data
        [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params); % Read data
        close all
        
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
        clear squid_cleaned

        % SQUID-MAG timelock
        params.modality = 'squid';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
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
    
    end

    %% Empty room & resting state for noise covariances
    if exist(fullfile(save_path, [params.sub '_resting_state_squid.mat']),'file') && overwrite.empty_room == false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    else
        clear data_ica
        data_ica = load(fullfile(save_path, [params.sub '_opm_timelocked.mat'])).timelocked{1};
        opm_chs = data_ica.label(contains(data_ica.label,'bz'));
        opm_grad = data_ica.grad;
        clear data_ica
        data_ica = load(fullfile(save_path, [params.sub '_squid_timelocked.mat'])).timelocked{1};
        squid_chs = data_ica.label(contains(data_ica.label,'MEG'));
        clear data_ica

        if any(strcmp(params.covs,'empty_room'))
            % Empty room
            opm_file = fullfile(raw_path, 'osmeg', 'EmptyRoomOPM_raw.fif');
            squid_file = fullfile(raw_path, 'meg', 'EmptyRoomMEG_tsss.fif');
            if exist(opm_file,'file') && exist(squid_file,'file')
                read_empty_rooms(opm_file, squid_file, opm_chs, squid_chs, opm_grad, save_path, params);
            end
        end
    
        if any(strcmp(params.covs,'resting_state'))
            % RESO
            opm_file = fullfile(raw_path, 'osmeg', 'RSEOOPM_raw.fif');
            aux_file = fullfile(raw_path, 'meg', 'RSEOEEG.fif');
            squid_file = fullfile(raw_path, 'meg', 'RSEOMEG_proc-tsss+corr98+mc+avgHead_meg.fif');
            if exist(opm_file,'file') && exist(squid_file,'file')
                read_osMEG_RS(opm_file, aux_file, opm_chs, save_path, params)
                read_cvMEG_RS(squid_file, squid_chs, save_path, params)
            end
        end

        clear squid_chs opm_grad opm_chs
    end

    %% HPI localization
    ft_hastoolbox('mne',1);
   
    if exist(fullfile(save_path_mri, 'opm_trans.mat'),'file') && overwrite.coreg==false
        disp(['Not overwriting OPM transform for ' params.sub]);
    else
        meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif']);
        if i_sub == 9
            meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_proc-tsss+corr98.fif']);
        end
        ft_hastoolbox('mne', 1);
        load(fullfile(save_path, [params.sub '_opm_ica_ds']));
        params.include_chs = data_ica_ds.label(find(contains(data_ica_ds.label,'bz')));
        opm_trans = fit_hpi(hpi_path, meg_file, save_path_mri, params);
        close all
        clear data_ica_ds

        %% Plot sensor arrays with headmodels and sourcemodel
        if exist(fullfile(save_path_mri, 'headmodels.mat'),'file') && exist(fullfile(save_path_mri, 'sourcemodel.mat'),'file') && ~isempty(opm_trans)
            opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelocked.mat'])).timelocked;
            opmeeg_timelockedT = load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat'])).timelocked;
            squideeg_timelocked = load(fullfile(save_path, [params.sub '_squideeg_timelocked.mat'])).timelocked;
            squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked.mat'])).timelocked;
            for i = 1:length(params.trigger_labels)
                opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.chanpos);
                opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.coilpos);
                opm_timelockedT{i}.grad.chanori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.chanori')';
                opm_timelockedT{i}.grad.coilori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.coilori')';
                opmeeg_timelockedT{i}.elec.chanpos = squideeg_timelocked{i}.elec.chanpos;
                opmeeg_timelockedT{i}.elec.elecpos = squideeg_timelocked{i}.elec.elecpos;
            end
            
            headmodel = load(fullfile(save_path_mri,'headmodels.mat')).headmodels.headmodel_meg;
            sourcemodel = load(fullfile(save_path_mri, 'sourcemodel')).sourcemodel;
            meshes = load(fullfile(save_path_mri,'meshes.mat')).meshes;  

            % Plot source and head models
            h=figure; 
            ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
            hold on; 
            ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
            ft_plot_headmodel(headmodel, 'facealpha', 0.25, 'edgealpha', 0.25)
            ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
            ft_plot_sens(opmeeg_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
            hold off;
            title('OPM-MEG')
            view([-140 10])
            savefig(h, fullfile(save_path, 'figs', 'opm_layout2.fig'))
            saveas(h, fullfile(save_path, 'figs', 'opm_layout2.jpg'))
            close all
    
            h=figure; 
            ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
            hold on; 
            ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.2,'FaceColor',[229 194 152]/256,'unit','cm')
            ft_plot_headmodel(headmodel, 'facealpha', 0.25, 'edgealpha', 0.25)
            ft_plot_sens(squid_timelocked{1}.grad,'unit','cm')
            ft_plot_sens(squideeg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
            hold off;
            title('SQUID-MEG')
            view([-140 10])
            savefig(h, fullfile(save_path, 'figs', 'meg_layout2.fig'))
            saveas(h, fullfile(save_path, 'figs', 'meg_layout2.jpg'))
            close all
            clear meshes sourcemodel headmodel squid_timelocked squideeg_timelocked
    
            %% Save
            save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
            save(fullfile(save_path, [params.sub '_opmeeg_timelockedT']), 'opmeeg_timelockedT', '-v7.3');
            clear opm_timelockedT opmeeg_timelockedT

        elseif isempty(opm_trans)
            warning(['HPI fit did not succeed for ' params.sub '.'])
            excl_subs_src = [excl_subs_src i_sub];
        else
            warning(['Headmodel/sourcemodel missing for ' params.sub '.'])
            excl_subs_src = [excl_subs_src i_sub];
        end
    end

    %% Dipole fits
    ft_hastoolbox('mne',1);  
    if exist(fullfile(save_path, [params.peaks{1}.label '_dipoles.mat']),'file') && overwrite.dip==false
        disp(['Not overwriting dipole source reconstruction for ' params.sub]);
    elseif exist(fullfile(save_path, [params.sub '_opm_timelockedT.mat']),'file')
        headmodel = load(fullfile(save_path_mri, 'headmodels.mat')).headmodels.headmodel_meg;
        mri_resliced = load(fullfile(save_path_mri, 'mri_resliced.mat')).mri_resliced;
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked.mat'])).timelocked;
        
        for i_peak = 1:length(params.peaks)
            peak_opm = load(fullfile(save_path, [params.sub '_opm_' params.peaks{i_peak}.label])).peak; 
            peak_squid = load(fullfile(save_path, [params.sub '_squid_' params.peaks{i_peak}.label])).peak; 
            fit_dipoles(save_path, squid_timelocked, opm_timelockedT, headmodel, mri_resliced, peak_squid, peak_opm, params);
            clear peak_opm peak_squid
        end
        clear squid_timelocked opm_timelockedT
    end

    %% Distributed source models
    ft_hastoolbox('mne',1);
    if exist(fullfile(save_path, 'opm_mne_peaks.mat'),'file') && overwrite.mne==false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    elseif exist(fullfile(save_path, [params.sub '_opm_timelockedT.mat']),'file')
        sourcemodel = load(fullfile(save_path_mri, 'sourcemodel')).sourcemodel;
        sourcemodel_inflated = load(fullfile(save_path_mri, 'sourcemodel_inflated')).sourcemodel_inflated;
        headmodel = load(fullfile(save_path_mri,'headmodels.mat')).headmodels.headmodel_meg;
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelockedT.mat'])).opm_timelockedT;
        squid_timelocked = load(fullfile(save_path, [params.sub '_squid_timelocked.mat'])).timelocked;
        
        for i = 1:length(params.trigger_codes)
            if any(strcmp(params.covs,'resting_state')) && exist(fullfile(save_path, [params.sub '_resting_state_squid.mat']),'file') && exist(fullfile(save_path, [params.sub '_resting_state_opm.mat']),'file')
                opm_timelockedT{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_opm.mat'])).opm_RS_cov;
                squid_timelocked{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_squid.mat'])).squid_RS_cov;
            end
            if any(strcmp(params.covs,'empty_room')) && exist(fullfile(save_path, [params.sub '_ER_squid.mat']),'file') && exist(fullfile(save_path, [params.sub '_ER_opm.mat']),'file')
                opm_timelockedT{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_opm.mat'])).opm_ER.cov;
                squid_timelocked{i}.cov_ER = load(fullfile(save_path, [params.sub '_ER_squid.mat'])).squid_ER.cov;
            end
        end
        
        %% MNE fit
        params.inv_method = 'mne';
        for i_cov = 1:length(params.covs)
            params.noise_cov = params.covs{i_cov}; 
            fit_mne(save_path, squid_timelocked, opm_timelockedT, headmodel, sourcemodel, sourcemodel_inflated, params);
        end

        %% Clear variables
        clear squid_timelocked opm_timelockedT headmodel sourcemodel sourcemodel_inflated
    end
end
close all

%% Save results in report
for i_sub = setdiff(subs_to_run,excl_subs)
    params.sub = ['sub_' num2str(i_sub,'%02d')];
    save_path = fullfile(base_save_path,params.paradigm,params.sub);
    create_sub_reports(save_path, i_sub, params);
end

%% Sensor level group analysis
if overwrite.sens_group
    save_path = fullfile(base_save_path,params.paradigm);
    if ~exist(fullfile(save_path,'figs'), 'dir')
           mkdir(fullfile(save_path,'figs'))
    end
    subs = setdiff(subs_to_run,excl_subs);
    sensor_results_goup(save_path,subs, params)
    close all
end

%% Dipole group analysis
if overwrite.dip_group
    save_path = fullfile(base_save_path,params.paradigm);
    subs = setdiff(subs_to_run,excl_subs_src);
    dipole_results_goup(save_path,subs, params)
end

%% MNE group analysis
if overwrite.mne_group
    save_path = fullfile(base_save_path,params.paradigm);
    subs = setdiff(subs_to_run,excl_subs_src);
    for i_cov = 1:length(params.covs)
        params.use_cov = params.covs{i_cov}; 
        mne_results_goup(save_path, subs, params);
    end
end

%% Group level report
if overwrite.mne_group
    save_path = fullfile(base_save_path,params.paradigm);
    subs = setdiff(subs_to_run,excl_subs_src);
    for i_peak = 1:length(params.peaks)
        create_group_report(save_path, subs, params, i_peak)
    end
end

%%
close all
clear all
% if on_server
%     exit
% end