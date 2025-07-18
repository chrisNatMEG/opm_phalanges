function [] = read_cvMEG_RS(squid_file, squid_chs, save_path, params)
%prprocess_cvMEG Read conventional MEG data for benchmarking
% recordings. 
% Requires the following arguments:
% Path: containing save_path and squid_file
% Params: containing pre, post (pre- and poststim).

%% --- Read triggers ---
cfg             = [];
cfg.datafile    = squid_file;
squid_raw         = ft_preprocessing(cfg);
squid_trig = contains(squid_raw.label,'STI101');
trig = squid_raw.trial{1}(squid_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trig = find(trig);
n_smpl = round((params.pre+params.post+2*params.pad)*squid_raw.fsample);
n_trl = floor((trig(end)-trig(1))/n_smpl -1);
trl_meg = zeros(n_trl,4);
trl_meg(:,1) = trig(1) + n_smpl*(0:(n_trl-1))';
trl_meg(:,2) = trig(1) + n_smpl*(1:n_trl)' - 1;
trl_meg(:,3) = -(params.pad+params.pre)*squid_raw.fsample;
trl_meg(:,4) = ones(length(trl_meg(:,1)),1);

% MEG data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
if ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
end
squid_RS_ica = ft_preprocessing(cfg, squid_raw);
clear squid_raw

cfg = [];
cfg.trl = trl_meg;
squid_RS_ica = ft_redefinetrial(cfg,squid_RS_ica);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
squid_RS_ica = ft_preprocessing(cfg,squid_RS_ica);

EOG_channels = find(contains(squid_RS_ica.label,'EOG'));
ECG_channels = find(contains(squid_RS_ica.label,'ECG'));
include_channels = [squid_chs; squid_RS_ica.label(EOG_channels); squid_RS_ica.label(ECG_channels)];

% Remove all but ExG and meg channel selection
cfg = [];
cfg.channel = include_channels;
squid_RS_ica = ft_selectdata(cfg,squid_RS_ica);

% MEG 
% Reject jump trials
cfg = [];
cfg.channel = 'meg';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,~] = ft_badsegment(cfg, squid_RS_ica);
squid_RS_ica = ft_rejectartifact(cfg,squid_RS_ica);

% Reject noisy trials
cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'std';
cfg.threshold = params.squidmag_std_threshold;
[cfg,~] = ft_badsegment(cfg, squid_RS_ica);
squid_RS_ica = ft_rejectartifact(cfg,squid_RS_ica);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'std';
cfg.threshold = params.squidgrad_std_threshold;
[cfg,~] = ft_badsegment(cfg, squid_RS_ica);
squid_RS_ica = ft_rejectartifact(cfg,squid_RS_ica);

cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'range';
cfg.threshold = params.squidmag_range_threshold;
[cfg,~] = ft_badsegment(cfg, squid_RS_ica);
squid_RS_ica = ft_rejectartifact(cfg,squid_RS_ica);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'range';
cfg.threshold = params.squidgrad_range_threshold;
[cfg,~] = ft_badsegment(cfg, squid_RS_ica);
squid_RS_ica = ft_rejectartifact(cfg,squid_RS_ica);

% Downsample
if isfield(params,'ds_freq') && ~isempty(params.ds_freq) && params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    squid_RS_ica = ft_resampledata(cfg, squid_RS_ica);
end


%% ICA
params.modality = 'squid';
params.layout = 'neuromag306mag.lay';
params.chs = 'meg';
squid_RS_ica = ica_MEG(squid_RS_ica, save_path, params, 0);

cfg = [];
cfg.channel = 'meg';
squid_RS_ica = ft_selectdata(cfg,squid_RS_ica);

%% Timelock
% Remove padding
cfg = [];
cfg.latency = [-params.pre params.post];
squid_RS_ica = ft_selectdata(cfg, squid_RS_ica);

% Demean
cfg = [];
cfg.demean = 'yes'; %% demean entire trial for whole trial cov
cfg.baselinewindow = [-params.pre 0];
squid_RS_ica = ft_preprocessing(cfg,squid_RS_ica);

% Average
cfg = [];
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'all';
squid_RS_cov = ft_timelockanalysis(cfg, squid_RS_ica).cov;

%% Save 
save(fullfile(save_path, [params.sub '_resting_state_squid']), 'squid_RS_ica', 'squid_RS_cov',"-v7.3");

end