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
    overwrite.preproc = false;
    overwrite.timelock = false;
    overwrite.TFR = false;
    overwrite.coreg = true;
    overwrite.mri = false;
    overwrite.dip = true;
    overwrite.empty_room = false;
    overwrite.mne = false;

    overwrite.sens_group = false;
    overwrite.dip_group = false;
    overwrite.mne_group = false;
else
    overwrite.preproc = true;
    overwrite.timelock = true;
    overwrite.TFR = true;
    overwrite.coreg = true;
    overwrite.mri = false;
    overwrite.dip = true;
    overwrite.empty_room = false;
    overwrite.mne = true;

    overwrite.sens_group = false;
    overwrite.dip_group = false;
    overwrite.mne_group = false;
end

%% Params
params = [];

params.paradigm = 'AudOdd';

% Trials
params.pre = 0.3; %0.1 sec
params.post = 0.75+0.2; %0.5 sec
params.pad = 0.2; %sec
params.delay = 0.01; % Stimulus delay in seconds (e.g., 0.01 for eartubes or 0.041 for membranes).
params.baseline = [-0.2 0];

% EEG
params.eeg_reref = 'all';%'EEG023';

% Filter
params.filter = [];
params.filter.hp_freq =1;%0.1;
params.filter.lp_freq = 20;%50;
params.filter.bp_freq = [];
params.filter.notch = [50 60]; %[50 60 100 120 150];

% Spatiotemporal filter (OPM-MEG only)
params.do_hfc = true;
params.hfc_order = 2;
params.do_amm = false;
params.amm_in = 12;
params.amm_out = 2;
params.amm_thr = 1;
params.do_ssp = false;
params.ssp_n = 6;

% Bad channel and trial detection thresholds
params.outlier_zscore = 3; % Outliers: how many stddevs above mean
params.outlier_ratio = 0.5; % Outliers: ratio of frequencies over threshold above which the channel is considered an outlier
params.corr_threshold = 0.6; % Correlation threshold for badchannei detection based on neighbor correlation
params.z_threshold = 20; % Zmax threshold for badchannel and trial detection based on jumps
params.opm_range_threshold = 20e-12; % Range for OPM badtrial detection
params.squidmag_range_threshold = 10e-12; % Range for SQUID-MAG badtrial detection
params.squidgrad_range_threshold = 4000e-13; % Range for SQUID-GRAD badtrial detection
params.debug = false; % Do manual rejection

% ICA ECG&EOG artifact removal 
params.n_comp = 40;
params.ica_cor = 0.8; % cutoff for EOG/ECG coherence
params.ica_coh = 0.95; % cutoff for EOG/ECG coherence

% Timelocking
params.ds_freq = 1000; % downsample frequency (timelock)
params.trigger_codes = {1 [3 11] [5 13]}; % combined oddball-nogo and oddball-go
params.trigger_labels = {'std' 'oddNoGo' 'oddGo'};
% params.trigger_codes = {1 3 5 11 13};
% params.trigger_labels = {'STD' 'LNG' 'LG' 'HNG' 'HG'};

% Evoked peaks to analyze
params.peaks = {};
params.peaks{1} = [];
params.peaks{1}.label = 'M100';
params.peaks{1}.peak_latency = [0.08 0.13];

% Time-Frequency
params.tfr = true;

% HPI coregistration
params.hpi_freq = 33;
params.hpi_gof = 0.9;

% Source reconstruction - dipoles
params.numdipoles = 2;

% Source reconstruction - distributed
params.source_fixedori = true; 
params.covs = {' '}; % noise cov to use; default=prestim, alt: 'resting_state', 'all', 'empty_room' , prestim = ' '
params.mne_view = 'sides';
params.plot_inflated = true;
params.target_region = {'superiortemporal', 'transversetemporal'};

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

bads =  {[]; %1
    {'R409_bz'}; %2
    {'L503_bz', 'R403_bz', 'R408_bz', 'R409_bz'}; %3 
    {'L109_bz', 'L604_bz', 'L205_bz', 'L209_bz', 'R403_bz', 'R408_bz'}; %4
    {'L102_bz', 'L310_bz', 'L209_bz', 'L505_bz', 'R403_bz', 'R408_bz', 'R409_bz'}; %5
    {'L409_bz', 'R304_bz', 'L310_bz', 'R403_bz', 'R408_bz', 'R409_bz'}; %6
    {'L209_bz', 'R403_bz', 'R408_bz', 'R409_bz'}; %7
    {'L209_bz', 'L502_bz', 'R209_bz', 'R403_bz', 'R402_bz', 'R502_bz', 'R409_bz'}; %8
    {'L502_bz', 'L503_bz', 'L209_bz', 'R408_bz'}; %9
    {'L111_bz', 'L209_bz', 'R408_bz'}; %10
    {'L603_bz', 'L209_bz', 'L109_bz', 'L111_bz', 'R408_bz'}; %11
    {'L604_bz', 'L209_bz', 'R408_bz'}; %12
    {'L604_bz', 'R401_bz', 'L209_bz', 'L109_bz', 'L111_bz', 'R408_bz'}; %13
    {'L209_bz', 'L411_bz', 'L504_bz', 'L401_bz', 'R409_bz'}; %14
    {'L502_bz', 'R403_bz', 'R409_bz'}}; %15

if on_server
    subs_to_run = 1:size(subses,1);
else
    subs_to_run = 4; %1:size(subses,1)
end
excl_subs = [3 15]; % split file TODO: allow split file
excl_subs_src = [1 excl_subs];

%% Loop over subjects
ssp_done = false(size(subs_to_run,2),1);
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
        if ~exist(meg_file,'file')
            meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_tsss.fif']);
        end
    end
    opm_file = fullfile(raw_path, 'osmeg', [params.paradigm 'OPM_raw.fif']);
    aux_file = fullfile(raw_path, 'meg', [params.paradigm 'EEG.fif']);
    params.ssp_file = fullfile(raw_path, 'osmeg', 'EmptyRoomOPM_raw.fif');
    
    %% Preprocessing
    if exist(fullfile(save_path, [params.sub '_opm_preproc.mat']),'file') && exist(fullfile(save_path, [params.sub '_squid_preproc.mat']),'file') && overwrite.preproc==false
        disp(['Not overwriting preproc for ' params.sub]);
    else
        ft_hastoolbox('mne', 1);

        %% OPM-MEG 
        % Read data
        [opm_cleaned, opmeeg_cleaned, ssp_done(i_sub)] = read_osMEG(opm_file, aux_file, save_path, params); % Read data
        close all

        % Correct old trigger codes
        if any(opm_cleaned.trialinfo==18)
            old_codes = [1 18 20 10 12];
            new_codes = {1 3 5 11 13};
            for i_code = 2:length(old_codes)
                opm_cleaned.trialinfo(opm_cleaned.trialinfo==old_codes(i_code)) = new_codes{i_code};
                opmeeg_cleaned.trialinfo(opmeeg_cleaned.trialinfo==old_codes(i_code)) = new_codes{i_code};
            end
        end

        % OPM ICA
        params.modality = 'opm';
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        data_ica = ica_MEG(opm_cleaned, save_path, params, 1);
        save(fullfile(save_path, [params.sub '_' params.modality '_preproc']), 'data_ica', '-v7.3'); 
        clear opm_cleaned data_ica

        % Make OPM-EEG layout
        cfg = [];
        cfg.elec = opmeeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_opmeeg_layout.mat']);
        opmeeg_layout = ft_prepare_layout(cfg);

        % OPM-EEG ICA
        params.layout = opmeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'opmeeg';
        data_ica = ica_MEG(opmeeg_cleaned, save_path, params, 1);
        save(fullfile(save_path, [params.sub '_' params.modality '_preproc']), 'data_ica', '-v7.3'); 
        clear opmeeg_cleaned data_ica
        
        %% SQUID-MEG 
        ft_hastoolbox('mne', 1);

        % Read data
        [squid_cleaned, squideeg_cleaned] = read_cvMEG(meg_file, save_path, params); % Read data
        
        % Correct old trigger codes
        if any(squid_cleaned.trialinfo==18)
            old_codes = [1 18 20 10 12];
            new_codes = {1 3 5 11 13};
            for i_code = 1:length(old_codes)
                squid_cleaned.trialinfo(squid_cleaned.trialinfo==old_codes(i_code)) = new_codes{i_code};
                squideeg_cleaned.trialinfo(squideeg_cleaned.trialinfo==old_codes(i_code)) = new_codes{i_code};
            end
        end

        % SQUID-MAG ICA
        params.modality = 'squid';
        params.layout = 'neuromag306mag.lay';
        params.chs = 'meg';
        data_ica = ica_MEG(squid_cleaned, save_path, params, 1); 
        save(fullfile(save_path, [params.sub '_' params.modality '_preproc']), 'data_ica', '-v7.3'); 
        clear squid_cleaned data_ica

        % Make SQUID-EEG layout
        cfg = [];
        cfg.elec = squideeg_cleaned.elec;
        cfg.output = fullfile(save_path, [params.sub '_squideeg_layout.mat']);
        megeeg_layout = ft_prepare_layout(cfg);

        % SQUID-EEG ICA
        params.layout = megeeg_layout;
        params.chs = 'EEG*';
        params.modality = 'squideeg';
        data_ica = ica_MEG(squideeg_cleaned, save_path, params, 1);
        save(fullfile(save_path, [params.sub '_' params.modality '_preproc']), 'data_ica', '-v7.3'); 
        clear squideeg_cleaned data_íca

        params = rmfield(params,{'modality', 'layout', 'chs'}); % remove fields used for picking modality
    end

    %% Timelocking
    if exist(fullfile(save_path, [params.sub '_opm_preproc.mat']),'file') && overwrite.timelock == false
        disp(['Not overwriting timelocked for ' params.sub]);
    else
        % OPM average
        params.modality = 'opm';
        opm_preproc = load(fullfile(save_path, [params.sub '_' params.modality '_preproc'])).data_ica; 
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        timelock_MEG(opm_preproc, save_path, params);
        clear opm_preproc

        % OPM-EEG Average
        params.modality = 'opmeeg';
        opmeeg_preproc = load(fullfile(save_path, [params.sub '_' params.modality '_preproc'])).data_ica; 
        params.layout = load(fullfile(save_path, [params.sub '_opmeeg_layout.mat'])).layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        timelock_MEG(opmeeg_preproc, save_path, params);
        clear opmeeg_preproc

        % SQUID-MAG timelock
        params.modality = 'squid';
        squid_preproc = load(fullfile(save_path, [params.sub '_' params.modality '_preproc'])).data_ica; 
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        timelock_MEG(squid_preproc, save_path, params);

        % SQUID-GRAD timelock
        squid_preproc = load(fullfile(save_path, [params.sub '_' params.modality '_preproc'])).data_ica; 
        params.modality = 'squidgrad';
        params.layout = 'neuromag306planar.lay';
        params.chs = 'megplanar';
        params.amp_scaler = 1;
        params.amp_label = 'B [T/cm]';
        timelock_MEG(squid_preproc, save_path, params);
        clear squid_preproc

        % SQUID-EEG timelock
        params.modality = 'squideeg';
        squideeg_preproc = load(fullfile(save_path, [params.sub '_' params.modality '_preproc'])).data_ica;
        params.layout = load(fullfile(save_path, [params.sub '_squideeg_layout.mat'])).layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        timelock_MEG(squideeg_preproc, save_path, params);
        clear squideeg_preproc

        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
    end

    %% TFR
    if exist(fullfile(save_path, [params.sub '_opm_tfr.mat']),'file') && overwrite.TFR == false
        disp(['Not overwriting TFR for ' params.sub]);
    else
        % OPM 
        params.modality = 'opm';
        opm_timelocked = load(fullfile(save_path, [params.sub '_opm_timelocked.mat'])).timelocked;
        params.layout = 'fieldlinebeta2bz_helmet.mat';
        params.chs = '*bz';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        TFR_MEG(opm_timelocked, save_path, params);

        % OPM-EEG
        params.modality = 'opmeeg';
        opmeeg_timelocked = load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat'])).timelocked;
        params.layout = load(fullfile(save_path, [params.sub '_opmeeg_layout.mat'])).layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        TFR_MEG(opmeeg_timelocked, save_path, params);

        % SQUID
        params.modality = 'squid';
        squid_timelocked = load(fullfile(save_path, [params.sub '_' params.modality '_timelocked'])).timelocked; 
        params.layout = 'neuromag306mag.lay';
        params.chs = 'megmag';
        params.amp_scaler = 1e15;
        params.amp_label = 'B [fT]';
        TFR_MEG(squid_timelocked, save_path, params);

        % SQUID-EEG
        params.modality = 'squideeg';
        squideeg_timelocked = load(fullfile(save_path, [params.sub '_' params.modality '_timelocked'])).timelocked;
        params.layout = load(fullfile(save_path, [params.sub '_squideeg_layout.mat'])).layout;
        params.chs = 'EEG*';
        params.amp_scaler = 1e9;
        params.amp_label = 'V [nV]';
        TFR_MEG(squideeg_timelocked, save_path, params);

        clear opm_timelocked opmeeg_timelocked squid_timelocked squideeg_timelocked
        params = rmfield(params,{'modality', 'layout', 'chs', 'amp_scaler', 'amp_label'}); % remove fields used for picking modality
    end

    %% Empty room & resting state for noise covariances
    if exist(fullfile(save_path, [params.sub '_resting_state_squid.mat']),'file') && overwrite.empty_room == false
        disp(['Not overwriting MNE source reconstruction for ' params.sub]);
    else
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
        opm_trans = load(fullfile(save_path_mri, 'opm_trans.mat')).opm_trans;
        opm_timelockedT = load(fullfile(save_path, [params.sub '_opm_timelocked.mat'])).timelocked;
        opmeeg_timelockedT = load(fullfile(save_path, [params.sub '_opmeeg_timelocked.mat'])).timelocked;
        for i = 1:length(params.trigger_labels)
            opm_timelockedT{i}.grad.chanpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.chanpos);
            opm_timelockedT{i}.grad.coilpos = opm_trans.transformPointsForward(opm_timelockedT{i}.grad.coilpos);
            opm_timelockedT{i}.grad.chanori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.chanori')';
            opm_timelockedT{i}.grad.coilori = (opm_trans.Rotation'*opm_timelockedT{i}.grad.coilori')';
            opmeeg_timelockedT{i}.elec.chanpos = squideeg_timelocked{i}.elec.chanpos;
            opmeeg_timelockedT{i}.elec.elecpos = squideeg_timelocked{i}.elec.elecpos;
        end
        save(fullfile(save_path, [params.sub '_opm_timelockedT']), 'opm_timelockedT', '-v7.3');
        save(fullfile(save_path, [params.sub '_opmeeg_timelockedT']), 'opmeeg_timelockedT', '-v7.3');
        clear opm_timelockedT opmeeg_timelockedT
    else
        meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_proc-tsss+corr98+mc+avgHead_meg.fif']);
        if ~exist(meg_file,'file')
            meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_proc-tsss+corr98.fif']);
            if ~exist(meg_file,'file')
                meg_file = fullfile(raw_path, 'meg', [params.paradigm 'MEG_tsss.fif']);
            end
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
            ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm')
            ft_plot_headmodel(headmodel, 'facealpha', 0.25, 'edgealpha', 0.25)
            ft_plot_sens(opm_timelockedT{1}.grad,'unit','cm')
            ft_plot_sens(opmeeg_timelockedT{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
            quiver3(0, 0, 0, 20, 0, 0,'k')
            text(20, 0, 0,'X')
            quiver3(0, 0, 0, 0, 20, 0,'k')
            text(0, 20, 0,'Y')
            quiver3(0, 0, 0, 0, 0, 20,'k')
            text(0, 0, 20,'Z')
            hold off;
            title('OPM-MEG')
            view([-140 10])
            savefig(h, fullfile(save_path, 'figs', 'opm_layout2.fig'))
            saveas(h, fullfile(save_path, 'figs', 'opm_layout2.jpg'))
            close all
    
            h=figure; 
            ft_plot_mesh(sourcemodel, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5,'unit','cm');
            hold on; 
            ft_plot_mesh(meshes(3),'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm')
            ft_plot_headmodel(headmodel, 'facealpha', 0.25, 'edgealpha', 0.25)
            ft_plot_sens(squid_timelocked{1}.grad,'unit','cm')
            ft_plot_sens(squideeg_timelocked{1}.elec,'unit','cm', 'style', '.r','elecsize',20)
            quiver3(0, 0, 0, 20, 0, 0,'k')
            text(20, 0, 0,'X')
            quiver3(0, 0, 0, 0, 20, 0,'k')
            text(0, 20, 0,'Y')
            quiver3(0, 0, 0, 0, 0, 20,'k')
            text(0, 0, 20,'Z')
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
        clear squid_timelocked opm_timelockedT mri_resliced headmodel
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
            if exist(fullfile(save_path, [params.sub '_resting_state_squid.mat']),'file') && exist(fullfile(save_path, [params.sub '_resting_state_opm.mat']),'file')
                opm_timelockedT{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_opm.mat'])).opm_RS_cov;
                squid_timelocked{i}.cov_RS = load(fullfile(save_path, [params.sub '_resting_state_squid.mat'])).squid_RS_cov;
            end
            if exist(fullfile(save_path, [params.sub '_ER_squid.mat']),'file') && exist(fullfile(save_path, [params.sub '_ER_opm.mat']),'file')
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
    params.do_sourcemovie = true;
    save_path = fullfile(base_save_path,params.paradigm);
    subs = setdiff(subs_to_run,excl_subs_src);
    save_path_mri = fullfile(base_save_path,'MRI',['sub_' num2str(subs(1),'%02d')]);
    sourcemodel = load(fullfile(save_path_mri, 'sourcemodel')).sourcemodel;
    sourcemodel_inflated = load(fullfile(save_path_mri, 'sourcemodel_inflated')).sourcemodel_inflated;
    for i_cov = 1:length(params.covs)
        params.use_cov = params.covs{i_cov}; 
        mne_results_goup(save_path, subs, sourcemodel, sourcemodel_inflated, params);
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
