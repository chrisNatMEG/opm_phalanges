function [] = read_osMEG_RS(opm_file, aux_file, opm_chs, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

%% --- Read triggers ---
% AUX
cfg = [];
cfg.datafile        = aux_file;
aux_raw = ft_preprocessing(cfg);
aux_trig = contains(aux_raw.label,'STI101');
trig = aux_raw.trial{1}(aux_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trig = find(trig);
aux_period = aux_raw.time{1}(trig(end))-aux_raw.time{1}(trig(1));
n_smpl = round((params.pre+params.post)*aux_raw.fsample);
n_trl = floor((trig(end)-trig(1))/n_smpl-1);
trl_aux = zeros(n_trl,4);
trl_aux(:,1) = trig(1) + n_smpl*(0:(n_trl-1))';
trl_aux(:,2) = trig(1) + n_smpl*(1:n_trl)' - 1;
trl_aux(:,3) = -params.pre*aux_raw.fsample;
trl_aux(:,4) = ones(length(trl_aux(:,1)),1);

% OPM
cfg = [];
cfg.datafile        = opm_file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
opm_raw = ft_preprocessing(cfg);
opm_trig = contains(opm_raw.label,'di');
trig = opm_raw.trial{1}(opm_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trig = find(trig);
opm_period = opm_raw.time{1}(trig(end))-opm_raw.time{1}(trig(1));
n_smpl = round((opm_period/aux_period)*(params.pre+params.post+2*params.pad)*opm_raw.fsample);
n_trl = floor((trig(end)-trig(1))/n_smpl-1);
trl_opm = zeros(n_trl,4);
trl_opm(:,1) = trig(1) + n_smpl*(0:(n_trl-1))';
trl_opm(:,2) = trig(1) + n_smpl*(1:n_trl)' - 1;
trl_opm(:,3) = -(params.pad+params.pre)*opm_raw.fsample;
trl_opm(:,4) = ones(length(trl_opm(:,1)),1);

% Check if uneven amount of trial. If so assume error in beginning.
if size(trl_aux,1) > size(trl_opm,1)
    trl_aux = trl_aux((end-size(trl_opm,1)+1):end,:);
elseif size(trl_aux,1) < size(trl_opm,1)
    trl_opm = trl_opm((end-size(trl_aux,1)+1):end,:);
end
if trl_aux(:,4) ~= trl_opm(:,4) % Throw error if trials don't match.
    disp()
    error('events do not match')
end

%% AUX data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
aux_epo = ft_preprocessing(cfg, aux_raw);

cfg = [];
cfg.trl = trl_aux;
aux_epo = ft_redefinetrial(cfg,aux_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
aux_epo = ft_preprocessing(cfg,aux_epo);

%% OPM data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
opm_epo = ft_preprocessing(cfg,opm_raw);

cfg = [];
cfg.trl             = trl_opm;
opm_epo = ft_redefinetrial(cfg,opm_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
opm_epo = ft_preprocessing(cfg,opm_epo);

%% --- Resample --- 
cfg            = [];
cfg.time = aux_epo.time;
cfg.detrend    = 'no';
opm_epo = ft_resampledata(cfg, opm_epo);

%% Select opm channels from main recording
cfg = [];
cfg.channel = opm_chs;
opm_epo = ft_selectdata(cfg, opm_epo);

%% Combine data
EOG_channels = find(contains(aux_epo.label,'EOG'));
ECG_channels = find(contains(aux_epo.label,'ECG'));
% EEG_channels = find(contains(aux_epo.label,'EEG'));
% MISC_channels = find(contains(aux_epo.label,'MISC'));
% TRIG_channels = find(contains(aux_epo.label,'STI101'));
include_channels = [EOG_channels; ECG_channels]; %; EEG_channels; MISC_channels; TRIG_channels];

comb = opm_epo; 
comb.elec = aux_epo.elec;
comb.time = aux_epo.time;
comb.label = [comb.label; aux_epo.label(include_channels)];
comb.hdr = aux_epo.hdr;
comb.hdr.label = comb.label;
comb.hdr.nChans = length(comb.label);
comb.hdr.chantype = [opm_epo.hdr.chantype; aux_epo.hdr.chantype(include_channels)];
comb.hdr.chanunit = [opm_epo.hdr.chanunit; aux_epo.hdr.chanunit(include_channels)];
comb.sampleinfo = aux_epo.sampleinfo;
comb.trialinfo = aux_epo.trialinfo;
n_smpl = size(aux_epo.trial{1},2);
for i = 1:length(comb.trial)
    comb.trial{i} = [comb.trial{i}(:,1:n_smpl); aux_epo.trial{i}(include_channels,:)]; 
end

%% Flip? 
chs = find(contains(comb.label,'bz'));
if isfield(params,'flip_sign') && params.flip_sign
    for i = 1:length(comb.trial)
        comb.trial{i}(chs,:) = -comb.trial{i}(chs,:);
    end
end

%% Reject bad channels
cfg = [];
cfg.trl = trl_opm;
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[~, ~, ~, ~, ~, ~, badtrl_opm_zmax] = opm_badchannels(cfg, opm_raw);

cfg = [];
cfg.trials  = setdiff(1:length(comb.trial),badtrl_opm_zmax); % remove bad trials
comb = ft_selectdata(cfg, comb);

cfg = []; % separate ExG channels
cfg.channel = {'EOG*', 'ECG*'};
ExG = ft_selectdata(cfg,comb);

%% Spatiotemporal filtering
if any(isnan(comb.grad.chanpos),'all')
    if isempty(setdiff(opm_chs,comb.grad.label))
        opm_grad.balance = comb.grad.balance;
        opm_grad.tra = eye(size(opm_grad.tra,1));
        comb.grad = opm_grad;
    else
        warning("OPM resting state: grad error. No covariance saved for %s",params.sub)
        return
    end
end

if params.do_hfc
    cfg = [];
    cfg.channel = '*bz';
    cfg.order = params.hfc_order;
    cfg.residualcheck = 'no';
    opm_cleaned = ft_denoise_hfc(cfg, comb);
elseif params.do_amm
    cfg = [];
    cfg.channel = '*bz';
    cfg.updatesens = 'yes';
    cfg.residualcheck = 'no';
    cfg.amm = [];
    cfg.amm.order_in = params.amm_in;
    cfg.amm.order_out = params.amm_out;
    cfg.amm.thr = params.amm_thr;
    opm_cleaned = ft_denoise_amm(cfg, comb);
else
    opm_cleaned = comb;
end

%% Recombine with ExG channels
opm_cleaned.label = vertcat(opm_cleaned.label,ExG.label);
opm_cleaned.hdr = comb.hdr;
incl = ismember(comb.hdr.label,opm_cleaned.label);
opm_cleaned.hdr.label = comb.hdr.label(incl);
opm_cleaned.hdr.chantype = comb.hdr.chantype(incl);
opm_cleaned.hdr.chanunit = comb.hdr.chanunit (incl);
for i = 1:length(opm_cleaned.trial)
    opm_cleaned.trial{i} = vertcat(opm_cleaned.trial{i}, ExG.trial{i}); 
end

%% Reject jump trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,~] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

%% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'std';
cfg.threshold = params.opm_std_threshold;
[cfg,~] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

%% Convert grad unit to cm to match TRIUX
opm_cleaned.grad = ft_convert_units(opm_cleaned.grad,'cm');

%% ICA
params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';
opm_RS_ica = ica_MEG(opm_cleaned, save_path, params, 0);

cfg = [];
cfg.channel = params.chs;
opm_RS_ica = ft_selectdata(cfg,opm_RS_ica);

%% Timelock
% Downsample
if params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    opm_RS_ica = ft_resampledata(cfg, opm_RS_ica);
end

% Remove padding
cfg = [];
cfg.latency = [-params.pre params.post];
opm_RS_ica = ft_selectdata(cfg, opm_RS_ica);

% Demean
cfg = [];
cfg.demean = 'yes'; %% demean entire trial for whole trial cov
opm_RS_ica = ft_preprocessing(cfg,opm_RS_ica);

% Average
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'all';
opm_RS_cov = ft_timelockanalysis(cfg, opm_RS_ica).cov;

%% Save
save(fullfile(save_path, [params.sub '_resting_state_opm']), 'opm_RS_ica', 'opm_RS_cov',"-v7.3");

end