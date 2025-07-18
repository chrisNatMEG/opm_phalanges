function read_empty_rooms(opm_file, squid_file, opm_chs, squid_chs, opm_grad, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

%% OPM
ft_hastoolbox('mne', 1);

% Read triggers
cfg = [];
cfg.datafile        = opm_file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
cfg.channel         = opm_chs;
data_raw = ft_preprocessing(cfg);
n_smpl = round((params.pre+params.post+2*params.pad)*data_raw.fsample);
n_trl = floor(data_raw.sampleinfo(2)/n_smpl-2);
if n_trl > 120
    n_trl = 120;
end
trl = zeros(n_trl,4);
trl(:,1) = data_raw.fsample + n_smpl*(0:(n_trl-1))';
trl(:,2) = data_raw.fsample + n_smpl*(1:n_trl)' - 1;
trl(:,3) = -(params.pad+params.pre)*data_raw.fsample;
trl(:,4) = ones(length(trl(:,1)),1);

% Filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
if ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
end
opm_ER_cleaned = ft_preprocessing(cfg, data_raw);

cfg = [];
cfg.trl = trl;
opm_ER_cleaned = ft_redefinetrial(cfg,opm_ER_cleaned);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
opm_ER_cleaned = ft_preprocessing(cfg,opm_ER_cleaned);

% Resample 
cfg                 = [];
cfg.method          = 'resample';
cfg.resamplefs      = 1000;
opm_ER_cleaned = ft_resampledata(cfg, opm_ER_cleaned);
opm_ER_cleaned.sampleinfo = round(trl(:,1:2)/5);

%% Reject bad channels
cfg = [];
cfg.trl = trl;
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[~, ~, ~, ~, ~, ~, badtrl_zmax] = opm_badchannels(cfg, data_raw);
clear data_raw

cfg = [];
cfg.trials  = setdiff(1:length(opm_ER_cleaned.trial),badtrl_zmax); % remove bad trials
opm_ER_cleaned = ft_selectdata(cfg, opm_ER_cleaned);

%% Spatiotemporal filtering
if any(isnan(opm_ER_cleaned.grad.chanpos),'all') 
    if isempty(setdiff(opm_chs,opm_ER_cleaned.grad.label))
        opm_grad.balance = opm_ER_cleaned.grad.balance;
        opm_grad.tra = eye(size(opm_grad.tra,1));
        opm_ER_cleaned.grad = opm_grad;
    else
        warning("OPM empty room: grad error. No covariance saved for %s",params.sub)
        return
    end
end

if params.do_hfc
    cfg = [];
    cfg.channel = '*bz';
    cfg.order = params.hfc_order;
    cfg.residualcheck = 'no';
    opm_ER_cleaned = ft_denoise_hfc(cfg, opm_ER_cleaned);
elseif params.do_amm
    cfg = [];
    cfg.channel = '*bz';
    cfg.updatesens = 'yes';
    cfg.residualcheck = 'no';
    cfg.amm = [];
    cfg.amm.order_in = params.amm_in;
    cfg.amm.order_out = params.amm_out;
    cfg.amm.thr = params.amm_thr;
    opm_ER_cleaned = ft_denoise_amm(cfg, opm_ER_cleaned);
else
    opm_ER_cleaned = opm_ER_cleaned;
end

%% Reject jump trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_jump] = ft_badsegment(cfg, opm_ER_cleaned);
opm_ER_cleaned = ft_rejectartifact(cfg,opm_ER_cleaned);

%% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'std';
cfg.threshold = params.opm_std_threshold;
[cfg,badtrl_std] = ft_badsegment(cfg, opm_ER_cleaned);
opm_ER_cleaned = ft_rejectartifact(cfg,opm_ER_cleaned);

cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'range';
cfg.threshold = params.opm_range_threshold;
[cfg,badtrl_range] = ft_badsegment(cfg, opm_ER_cleaned);
opm_ER_cleaned = ft_rejectartifact(cfg,opm_ER_cleaned);

[~,idx]=ismember(opm_ER_cleaned.sampleinfo,badtrl_jump,'rows');
badtrl_jump = find(idx);
[~,idx]=ismember(opm_ER_cleaned.sampleinfo,badtrl_std,'rows');
badtrl_std = find(idx);
[~,idx]=ismember(opm_ER_cleaned.sampleinfo,badtrl_range,'rows');
badtrl_range = find(idx);

save(fullfile(save_path, [params.sub '_opm_badtrls']), ...
    'badtrl_jump', ...
    'badtrl_std', ...
    'badtrl_range', ...
    'badtrl_zmax' ,"-v7.3"); 

% Flip?
chs = find(contains(opm_ER_cleaned.label,'bz'));
if isfield(params,'flip_sign') && params.flip_sign
    for i = 1:length(opm_ER_cleaned.trial)
        opm_ER_cleaned.trial{i}(chs,:) = -opm_ER_cleaned.trial{i}(chs,:);
    end
end

%% Downsample
if isfield(params,'ds_freq') && ~isempty(params.ds_freq) && params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    opm_ER_cleaned = ft_resampledata(cfg, ppm_ER_cleaned);
end

%% Timelock

% Remove padding
cfg = [];
cfg.channel = '*bz';
cfg.latency = [-params.pre params.post];
opm_ER_cleaned = ft_selectdata(cfg, opm_ER_cleaned);

% Demean
cfg = [];
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
opm_ER_cleaned = ft_preprocessing(cfg, opm_ER_cleaned);

% Average
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'all';
opm_ER_tlk = ft_timelockanalysis(cfg, opm_ER_cleaned);

%% Save
save(fullfile(save_path, [params.sub '_ER_opm']), 'opm_ER_cleaned', 'opm_ER_tlk',"-v7.3");
clear opm_ER_cleaned opm_ER_cleaned

%% SQUID
ft_hastoolbox('mne', 1);

% Read triggers
cfg             = [];
cfg.datafile    = squid_file;
cfg.channel         = squid_chs;
cfg.checkmaxfilter = 'no';
data_raw         = ft_preprocessing(cfg);
n_smpl = round((params.pre+params.post+2*params.pad)*data_raw.fsample);
n_trl = floor(data_raw.sampleinfo(2)/n_smpl -2);
if n_trl > 120
    n_trl = 120;
end
trl = zeros(n_trl,4);
trl(:,1) = data_raw.fsample + n_smpl*(0:(n_trl-1))';
trl(:,2) = data_raw.fsample + n_smpl*(1:n_trl)' - 1;
trl(:,3) = -(params.pad+params.pre)*data_raw.fsample;
trl(:,4) = ones(length(trl(:,1)),1);

% Filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
if ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
end
squid_ER_cleaned = ft_preprocessing(cfg, data_raw);
clear data_raw

cfg = [];
cfg.trl = trl;
squid_ER_cleaned = ft_redefinetrial(cfg,squid_ER_cleaned);
clear trl 

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
squid_ER_cleaned = ft_preprocessing(cfg,squid_ER_cleaned);

% Reject jump trials
cfg = [];
cfg.channel = 'meg';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_squid_jump] = ft_badsegment(cfg, opm_ER_cleaned);
squid_ER_cleaned = ft_rejectartifact(cfg,opm_ER_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'std';
cfg.threshold = params.squidmag_std_threshold;
[cfg,badtrl_squidmag_std] = ft_badsegment(cfg, squid_ER_cleaned);
squid_ER_cleaned = ft_rejectartifact(cfg,squid_ER_cleaned);

cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'range';
cfg.threshold = params.squidmag_range_threshold;
[cfg,badtrl_squidmag_range] = ft_badsegment(cfg, squid_ER_cleaned);
squid_ER_cleaned = ft_rejectartifact(cfg,squid_ER_cleaned);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'std';
cfg.threshold = params.squidgrad_std_threshold;
[cfg,badtrl_squidgrad_std] = ft_badsegment(cfg, squid_ER_cleaned);
squid_ER_cleaned = ft_rejectartifact(cfg,squid_ER_cleaned);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'range';
cfg.threshold = params.squidgrad_range_threshold;
[cfg,badtrl_squidgrad_range] = ft_badsegment(cfg, squid_ER_cleaned);
squid_ER_cleaned = ft_rejectartifact(cfg,squid_ER_cleaned);

[~,idx]=ismember(squid_ER_cleaned.sampleinfo,badtrl_squid_jump,'rows');
badtrl_squid_jump = find(idx);
[~,idx]=ismember(squid_ER_cleaned.sampleinfo,badtrl_squidmag_std,'rows');
badtrl_squid_std = find(idx);
[~,idx]=ismember(squid_ER_cleaned.sampleinfo,badtrl_squidgrad_std,'rows');
badtrl_squid_std = unique([badtrl_squid_std; find(idx)]);
[~,idx]=ismember(squid_ER_cleaned.sampleinfo,badtrl_squidmag_range,'rows');
badtrl_squid_range = find(idx);
[~,idx]=ismember(squid_ER_cleaned.sampleinfo,badtrl_squidgrad_range,'rows');
badtrl_squid_range = unique([badtrl_squid_range; find(idx)]);
save(fullfile(save_path, [params.sub '_squid_badtrls']), ...
    'badtrl_squid_jump', ...
    'badtrl_squid_range', ...
    'badtrl_squid_std',"-v7.3"); 

%% Downsample
if isfield(params,'ds_freq') && ~isempty(params.ds_freq) && params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    squid_ER_cleaned = ft_resampledata(cfg, squid_ER_cleaned);
end

%% Timelock
% Remove padding
cfg = [];
cfg.channel = 'meg';
cfg.latency = [-params.pre params.post];
squid_ER_cleaned = ft_selectdata(cfg, squid_ER_cleaned);

% Demean
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-params.pre 0];
squid_ER_cleaned = ft_preprocessing(cfg, squid_ER_cleaned);

% Average
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'all';
squid_ER_cov = ft_timelockanalysis(cfg, squid_ER_cleaned).cov;

%% Save
save(fullfile(save_path, [params.sub '_ER_squid']), 'squid_ER_cleaned', 'squid_ER_cov', "-v7.3");

clear opm_ER_cleaned data_raw opm_ER_cleaned
end