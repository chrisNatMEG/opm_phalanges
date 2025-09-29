function [squid_cleaned, squideeg_cleaned] = read_cvMEG(squid_file, save_path, params)
%prprocess_cvMEG Read conventional MEG data for benchmarking
% recordings. 
% Requires the following arguments:
% Path: containing save_path and squid_file
% Params: containing pre, post (pre- and poststim).

%% --- Read triggers ---
trl_squid = [];
cfg             = [];
cfg.datafile    = squid_file;
squid_raw         = ft_preprocessing(cfg);
squid_trig = find(contains(squid_raw.label,'STI101'));
trig = squid_raw.trial{1}(squid_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl_squid(:,1) = find(trig)-(params.pre+params.pad)*squid_raw.fsample;
trl_squid(:,2) = find(trig)+(params.post+params.pad)*squid_raw.fsample;
trl_squid(:,3) = -(params.pre+params.pad)*squid_raw.fsample;
trl_squid(:,4) = squid_raw.trial{1}(squid_trig,trig);
trl_squid(:,1:2) = trl_squid(:,1:2) + floor(params.delay*squid_raw.fsample); % adjust for stim delay
trl_squid = round(trl_squid);

%% MEG data filter & epoch
cfg = [];
if isfield(params.filter,'lp_freq') && ~isempty(params.filter.lp_freq)
    cfg.lpfilter        = 'yes';         
    cfg.lpfreq          = params.filter.lp_freq;
end
if isfield(params.filter,'hp_freq') && ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
end
if isfield(params.filter,'bp_freq') && ~isempty(params.filter.bp_freq)
    cfg.bpfilter        = 'yes'; 
    cfg.bpfreq          = params.filter.bp_freq;
    cfg.bpinstabilityfix  = 'reduce';
end
squid_epo = ft_preprocessing(cfg,squid_raw);

cfg = [];
cfg.trl             = trl_squid;
squid_epo = ft_redefinetrial(cfg,squid_epo);

cfg = [];
cfg.dftfilter    = 'yes';        
cfg.dftfreq      = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
squid_epo = ft_preprocessing(cfg,squid_epo);

% Find bad channels
cfg = [];
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[badchs_opmeeg, badchs_squideeg_flat, badchs_squideeg_neighbors] = eeg_badchannels(cfg,squid_raw);
clear squid_raw

%% MEG 
cfg = [];
cfg.channel = squid_epo.label(find(~contains(squid_epo.label,'eeg')));
squid_cleaned = ft_selectdata(cfg, squid_epo);

% no bad channel detection since maxfilter already does that

% Reject jump trials
cfg = [];
cfg.channel = 'meg';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_squid_jump] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

badtrl_squidmag_std = zeros(0,2);
badtrl_squidgrad_std = zeros(0,2);
% Reject noisy trials
% cfg = [];
% cfg.channel = 'megmag';
% cfg.metric = 'std';
% cfg.threshold = params.squidmag_std_threshold;
% [cfg,badtrl_squidmag_std] = ft_badsegment(cfg, squid_cleaned);
% squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);
% 
% cfg = [];
% cfg.channel = 'megplanar';
% cfg.metric = 'std';
% cfg.threshold = params.squidgrad_std_threshold;
% [cfg,badtrl_squidgrad_std] = ft_badsegment(cfg, squid_cleaned);
% squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

cfg = [];
cfg.channel = 'megmag';
cfg.metric = 'range';
cfg.threshold = params.squidmag_range_threshold;
[cfg,badtrl_squidmag_range] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

cfg = [];
cfg.channel = 'megplanar';
cfg.metric = 'range';
cfg.threshold = params.squidgrad_range_threshold;
[cfg,badtrl_squidgrad_range] = ft_badsegment(cfg, squid_cleaned);
squid_cleaned = ft_rejectartifact(cfg,squid_cleaned);

%% EEG
cfg = [];
cfg.channel = squid_epo.label(find(~contains(squid_epo.label,'MEG')));
squideeg_cleaned = ft_selectdata(cfg, squid_epo);
clear squid_epo

cfg = [];
cfg.channel = {'EOG', 'ECG'};
exg = ft_selectdata(cfg, squideeg_cleaned);

% Interpolate bad chs
cfg = [];
cfg.method = 'triangulation';
cfg.senstype = 'EEG';
neighbors = ft_prepare_neighbours(cfg,squideeg_cleaned);
cfg = [];
cfg.method = 'spline';
cfg.neighbors = neighbors;
cfg.badchannel = badchs_opmeeg;
cfg.senstype = 'EEG';
squideeg_cleaned = ft_channelrepair(cfg, squideeg_cleaned);

% Re-reference
cfg = [];
cfg.channel = 'all';
cfg.refef = 'yes';
cfg.refmethod = 'avg';
cfg.reffchannel = 'all';
squideeg_cleaned = ft_preprocessing(cfg,squideeg_cleaned);

% Append ECG/EOG channels
cfg = [];
squideeg_cleaned = ft_appenddata(cfg,squideeg_cleaned,exg);

% Reject jump trials
cfg = [];
cfg.channel = {'EEG*'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg, badtrl_squideeg_jump] = ft_badsegment(cfg, squideeg_cleaned);
squideeg_cleaned = ft_rejectartifact(cfg,squideeg_cleaned);

badtrl_squideeg_std = zeros(0,2);
% % Reject noisy trials
% cfg = [];
% cfg.channel = {'EEG*'};
% cfg.metric = 'std';
% cfg.threshold = params.eeg_std_threshold;
% [cfg, badtrl_squideeg_std] = ft_badsegment(cfg, squideeg_cleaned);
% squideeg_cleaned = ft_rejectartifact(cfg,squideeg_cleaned);

%% Save 
save(fullfile(save_path, [params.sub '_squideeg_badchs']), ...
    'badchs_squideeg_flat', ...
    'badchs_squideeg_neighbors',"-v7.3"); 

[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squid_jump,'rows');
badtrl_squid_jump = find(idx);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidmag_std,'rows');
badtrl_squid_std = find(idx);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidgrad_std,'rows');
badtrl_squid_std = unique([badtrl_squid_std; find(idx)]);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidmag_range,'rows');
badtrl_squid_range = find(idx);
[~,idx]=ismember(squid_cleaned.sampleinfo,badtrl_squidgrad_range,'rows');
badtrl_squid_range = unique([badtrl_squid_std; find(idx)]);
save(fullfile(save_path, [params.sub '_squid_badtrls']), ...
    'badtrl_squid_jump', ...
    'badtrl_squid_jump', ...
    'badtrl_squid_std',"-v7.3"); 

[~,idx]=ismember(squideeg_cleaned.sampleinfo,badtrl_squideeg_jump,'rows');
badtrl_squideeg_jump = find(idx);
[~,idx]=ismember(squideeg_cleaned.sampleinfo,badtrl_squideeg_std,'rows');
badtrl_squideeg_std = find(idx);
save(fullfile(save_path, [params.sub '_squideeg_badtrls']), ...
    'badtrl_squideeg_jump', ...
    'badtrl_squideeg_std', "-v7.3"); 

%save(fullfile(save_path, [params.sub '_squid_cleaned']), 'squid_cleaned',"-v7.3");
%save(fullfile(save_path, [params.sub '_squideeg_cleaned']), 'squideeg_cleaned',"-v7.3"); disp('done');

%% Downsample
if isfield(params,'ds_freq') && ~isempty(params.ds_freq) && params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    squid_cleaned = ft_resampledata(cfg, squid_cleaned);
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    squideeg_cleaned = ft_resampledata(cfg, squideeg_cleaned);
end

%% Remove padding
% cfg = [];
% cfg.latency = [-params.pre params.post];
% squid_cleaned = ft_selectdata(cfg, squid_cleaned); 
% 
% cfg = [];
% cfg.latency = [-params.pre params.post];
% squideeg_cleaned = ft_selectdata(cfg, squideeg_cleaned); 

end