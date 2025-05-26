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
n_smpl = round((params.pre+params.post)*squid_raw.fsample);
n_trl = floor((trig(end)-trig(1))/n_smpl -1);
trl_meg = zeros(n_trl,4);
trl_meg(:,1) = trig(1) + n_smpl*(0:(n_trl-1))';
trl_meg(:,2) = trig(1) + n_smpl*(1:n_trl)' - 1;
trl_meg(:,3) = -params.pre*squid_raw.fsample;
trl_meg(:,4) = ones(length(trl_meg(:,1)),1);

%% MEG data filter & epoch
cfg = [];
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
cfg.hpfilter        = 'yes';         
cfg.hpfreq          = params.filter.hp_freq;
cfg.hpinstabilityfix  = 'reduce';
squid_epo = ft_preprocessing(cfg, squid_raw);

cfg = [];
cfg.trl = trl_meg;
squid_epo = ft_redefinetrial(cfg,squid_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
squid_epo = ft_preprocessing(cfg,squid_epo);

EOG_channels = find(contains(squid_epo.label,'EOG'));
ECG_channels = find(contains(squid_epo.label,'ECG'));
include_channels = [squid_chs; squid_epo.label(EOG_channels); squid_epo.label(ECG_channels)];

% Remove all but ExG and meg channel selection
cfg = [];
cfg.channel = include_channels;
squid_epo = ft_selectdata(cfg,squid_epo);

%% MEG 
% Reject jump trials
cfg = [];
cfg.channel = 'meg';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,~] = ft_badsegment(cfg, squid_epo);
squid_cleaned = ft_rejectartifact(cfg,squid_epo);

% Reject noisy trials
cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'std';
cfg.threshold = params.squidmag_std_threshold;
[cfg,~] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'std';
cfg.threshold = params.squidgrad_std_threshold;
[cfg,~] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

%% ICA
params.modality = 'squid';
params.layout = 'neuromag306mag.lay';
params.chs = 'meg';
squid_RS_ica = ica_MEG(squid_cleaned, save_path, params, 0);

cfg = [];
cfg.channel = 'meg';
squid_RS_ica = ft_selectdata(cfg,squid_RS_ica);

%% Timelock
% Downsample
if params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    squid_RS_ica = ft_resampledata(cfg, squid_RS_ica);
end

% Remove padding
cfg = [];
cfg.channel = params.chs;
cfg.latency = [-params.pre params.post];
squid_RS_ica = ft_selectdata(cfg, squid_RS_ica);

% Demean
cfg = [];
cfg.demean = 'yes'; %% demean entire trial for whole trial cov
squid_RS_ica = ft_preprocessing(cfg,squid_RS_ica);

% Average
cfg = [];
cfg.channel = 'megmag';
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'all';
squidmag_RS_tlk = ft_timelockanalysis(cfg, squid_RS_ica);
cfg = [];
cfg.channel = 'megplanar';
cfg.covariance          = 'yes';
cfg.covariancewindow    = 'all';
squidgrad_RS_tlk = ft_timelockanalysis(cfg, squid_RS_ica);

%% Save 
save(fullfile(save_path, [params.sub '_resting_state_squid']), 'squid_RS_ica', 'squidmag_RS_tlk', 'squidgrad_RS_tlk',"-v7.3");

end