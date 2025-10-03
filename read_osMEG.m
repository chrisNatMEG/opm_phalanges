function [opm_cleaned, opmeeg_cleaned, ssp_done] = read_osMEG(opm_file, aux_file, save_path, params)
%prprocess_osMEG Read on-scalp MEG data for benchmarking
% recordings and combine with auxiliary TRIUX data/EEG. 
% Requires the following arguments:
% Path: containing save_path and meg_file
% Params: containing pre, post (pre- and poststim), and ds_freq 
% (downsampling frequency).

%% --- Read triggers ---
% OPM
trl_opm=[];
cfg = [];
cfg.datafile        = opm_file;
cfg.coordsys        = 'dewar';
cfg.coilaccuracy    = 0;
opm_raw = ft_preprocessing(cfg);
opm_trig = find(contains(opm_raw.label,'di'));
trig = opm_raw.trial{1}(opm_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl_opm(:,1) = find(trig)-(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,2) = find(trig)+(params.post+params.pad)*opm_raw.fsample;
trl_opm(:,3) = -(params.pre+params.pad)*opm_raw.fsample;
trl_opm(:,4) = opm_raw.trial{1}(opm_trig,trig);
trl_opm(:,1:2) = trl_opm(:,1:2) + floor(params.delay*opm_raw.fsample); % adjust for stim delay
trl_opm = round(trl_opm);

% AUX
trl_aux=[];
cfg = [];
cfg.datafile        = aux_file;
aux_raw = ft_preprocessing(cfg);
aux_trig = find(contains(aux_raw.label,'STI101'));
trig = aux_raw.trial{1}(aux_trig,:)>0.5;
trig = [false trig(2:end)&~trig(1:end-1)];
trl_aux(:,1) = find(trig)-(params.pre+params.pad)*aux_raw.fsample;
trl_aux(:,2) = find(trig)+(params.post+params.pad)*aux_raw.fsample;
trl_aux(:,3) = -(params.pre+params.pad)*aux_raw.fsample;
trl_aux(:,4) = aux_raw.trial{1}(aux_trig,trig);
trl_aux(:,1:2) = trl_aux(:,1:2) + floor(params.delay*aux_raw.fsample); % adjust for stim delay
trl_aux = round(trl_aux);

% Check if uneven amount of trial. If so assume error in beginning.
if size(trl_aux,1) > size(trl_opm,1)
    trl_aux = trl_aux((end-size(trl_opm,1)+1):end,:);
elseif size(trl_aux,1) < size(trl_opm,1)
    trl_opm = trl_opm((end-size(trl_aux,1)+1):end,:);
end
if trl_aux(:,4) ~= trl_opm(:,4) % Throw error if trials don't match.
    error('events do not match')
end

%% AUX data filter & epoch
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
aux_epo = ft_preprocessing(cfg,aux_raw);

cfg = [];
cfg.trl             = trl_aux;
aux_epo = ft_redefinetrial(cfg,aux_epo);

cfg = [];
cfg.dftfilter       = 'yes';        
cfg.dftfreq         = params.filter.notch;
aux_epo = ft_preprocessing(cfg,aux_epo);

% Find bad channels
cfg = [];
cfg.trl = trl_aux;
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[badchs_opmeeg, badchs_opmeeg_flat, badchs_opmeeg_neighbors] = eeg_badchannels(cfg,aux_raw);
clear aux_raw

%% OPM data filter & epoch
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
opm_epo = ft_preprocessing(cfg, opm_raw);

cfg = [];
cfg.trl             = trl_opm;
opm_epo = ft_redefinetrial(cfg,opm_epo);

cfg = [];
cfg.dftfilter       = 'yes';        
cfg.dftfreq         = params.filter.notch;
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
opm_epo = ft_preprocessing(cfg, opm_epo);

% Find bad opm channels
[badchs_opm, badchs_opm_flat, badchs_opm_std, badchs_opm_neighbors, badchs_opm_outlier, badchs_opm_noloc] = opm_badchannels(opm_raw,trl_opm,params,save_path);
clear opm_raw

%% --- Resample --- 
cfg            = [];
cfg.time = aux_epo.time;
cfg.detrend    = 'no';
opm_epo = ft_resampledata(cfg, opm_epo);

%% Combine data
EOG_channels = find(contains(aux_epo.label,'EOG'));
ECG_channels = find(contains(aux_epo.label,'ECG'));
EEG_channels = find(contains(aux_epo.label,'EEG'));
MISC_channels = find(contains(aux_epo.label,'MISC'));
TRIG_channels = find(contains(aux_epo.label,'STI101'));
include_channels = [EOG_channels; ECG_channels; EEG_channels; MISC_channels; TRIG_channels];

comb = opm_epo; 
comb.elec = aux_epo.elec;
comb.time = aux_epo.time;
comb.label = [opm_epo.label; aux_epo.label(include_channels)];
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
clear opm_epo aux_epo

%% Flip? 
chs = find(contains(comb.label,'bz'));
if isfield(params,'flip_sign') && params.flip_sign
    for i = 1:length(comb.trial)
        comb.trial{i}(chs,:) = -comb.trial{i}(chs,:);
    end
    comb.grad.chanori = -comb.grad.chanori;
    comb.grad.coilori = -comb.grad.coilori;
end

%% OPM 
if ~isempty(params.manual_bads)
    badchs_opm = [badchs_opm; params.manual_bads];
    badchs_opm_manual = params.manual_bads;
else
    badchs_opm_manual = [];
end

cfg = [];
cfg.channel = setdiff(comb.label,badchs_opm);
comb = ft_selectdata(cfg, comb);

cfg = []; % separate ExG channels
cfg.channel = {'EOG*', 'ECG*'};
ExG = ft_selectdata(cfg,comb);

%% Spatiotemporal filtering
if params.do_ssp && isfile(params.ssp_file)
    ssp_done = true;
    % Load reference data
    cfg = [];
    cfg.datafile        = params.ssp_file;
    cfg.coordsys        = 'dewar';
    cfg.coilaccuracy    = 0;
    refdata = ft_preprocessing(cfg);
    % Create trials (same length as main data)
    n_smpl = round((params.pre+params.post+2*params.pad)*refdata.fsample);
    n_trl = floor(refdata.sampleinfo(2)/n_smpl-2);
    trl_ref = zeros(n_trl,4);
    trl_ref(:,1) = refdata.fsample + n_smpl*(0:(n_trl-1))';
    trl_ref(:,2) = refdata.fsample + n_smpl*(1:n_trl)' - 1;
    trl_ref(:,3) = -(params.pad+params.pre)*refdata.fsample;
    trl_ref(:,4) = ones(length(trl_ref(:,1)),1);
    % Filter data (same as OPM data)
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
    refdata = ft_preprocessing(cfg, refdata);
    % Segment
    cfg = [];
    cfg.trl             = trl_ref;
    refdata = ft_redefinetrial(cfg,refdata);
    % Notch filter
    cfg = [];
    cfg.dftfilter       = 'yes';        
    cfg.dftfreq         = params.filter.notch;
    refdata = ft_preprocessing(cfg, refdata);
    % Downsample
    cfg            = [];
    cfg.resamplefs = comb.fsample;
    cfg.detrend    = 'no';
    refdata = ft_resampledata(cfg, refdata);
    % Select OPMs present in main data
    cfg = [];
    cfg.channel         = '*_b*';
    opm_cleaned = ft_selectdata(cfg,comb);
    datachans = opm_cleaned.label;
    cfg = [];
    cfg.channel         = datachans;
    refdata = ft_selectdata(cfg,refdata);
    refchans = refdata.label;
    if ~isempty(setdiff(datachans,refchans))
        ft_warning('SSP error: ref is missing channels present in data.');
        ssp_done = false;
    end
    [~, i_refchans, i_datachans] = intersect(refchans,datachans,'stable');
    % Calculate PCs and construct projector
    [coeff, ~, ~, ~, ~] = pca(cell2mat(refdata.trial)','NumComponents',params.ssp_n);
    proj = eye(size(coeff,1)) - coeff*transpose(coeff);
    % Apply projector    
    for i_trl = 1:length(opm_cleaned.trial)
        opm_cleaned.trial{i_trl}(i_datachans,:) = proj * opm_cleaned.trial{i_trl}(i_datachans,:);
    end
else
    ssp_done = false;
end
if params.do_hfc && ~ssp_done
    cfg = [];
    cfg.channel = '*_b*';
    cfg.order = params.hfc_order;
    cfg.residualcheck = 'no';
    opm_cleaned = ft_denoise_hfc(cfg, comb);
elseif params.do_amm && ~ssp_done
    cfg = [];
    cfg.channel = '*_b*';
    cfg.updatesens = 'yes';
    cfg.residualcheck = 'no';
    cfg.amm = [];
    cfg.amm.order_in = params.amm_in;
    cfg.amm.order_out = params.amm_out;
    cfg.amm.thr = params.amm_thr;
    opm_cleaned = ft_denoise_amm(cfg, comb);
elseif ssp_done
    disp('SSP successuful');
else
    cfg = [];
    cfg.channel         = '*_b*';
    opm_cleaned = ft_selectdata(cfg,comb);
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

% Reject jump trials
smplinfo = opm_cleaned.sampleinfo;

if isfield(params,'debug') && params.debug == 1
    cfg = [];
    cfg.method = 'summary';
    cfg.channel = '*bz';
    dummy = ft_rejectvisual(cfg,opm_cleaned);
end

cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_opm_jump] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

badtrl_opm_std = zeros(0,2);
% % Reject noisy trials
% cfg = [];
% cfg.channel = {'*bz'};
% cfg.metric = 'std';
% cfg.threshold = params.opm_std_threshold;
% [cfg,badtrl_opm_std] = ft_badsegment(cfg, opm_cleaned);
% opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'range';
cfg.threshold = params.opm_range_threshold;
[cfg,badtrl_opm_range] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

% Convert grad unit to cm to match TRIUX grad
opm_cleaned.grad = ft_convert_units(opm_cleaned.grad,'cm');

% Save bad channels 
save(fullfile(save_path, [params.sub '_opm_badchs']), ...
    'badchs_opm_flat', ...
    'badchs_opm_std', ...
    'badchs_opm_neighbors', ...
    'badchs_opm_outlier' , ...
    'badchs_opm_manual' , ...
    'badchs_opm_noloc',"-v7.3"); 

% Save bad trials
[~,idx]=ismember(smplinfo,badtrl_opm_jump,'rows');
badtrl_opm_jump = find(idx);
[~,idx]=ismember(smplinfo,badtrl_opm_std,'rows');
badtrl_opm_std = find(idx);
[~,idx]=ismember(smplinfo,badtrl_opm_range,'rows');
badtrl_opm_range = find(idx);
save(fullfile(save_path, [params.sub '_opm_badtrls']), ...
    'badtrl_opm_jump', ...
    'badtrl_opm_std', ...
    'badtrl_opm_range',"-v7.3"); 

%% EEG
cfg = [];
cfg.channel = 'EEG';
opmeeg_cleaned = ft_selectdata(cfg, comb);
clear comb

% Interpolate bad chs
% cfg = [];
% cfg.method = 'triangulation';
% cfg.senstype = 'EEG';
% neighbors = ft_prepare_neighbours(cfg,opmeeg_cleaned);
% cfg = [];
% cfg.method = 'spline';
% cfg.neighbors = neighbors;
% cfg.badchannel = badchs_opmeeg;
% cfg.senstype = 'EEG';
% opmeeg_cleaned = ft_channelrepair(cfg, opmeeg_cleaned);

% Re-reference
if isfield(params,'eeg_reref') && strcmp(params.eeg_reref,'avg')
    cfg = [];
    cfg.channel = 'all';
    cfg.reref = 'yes';
    cfg.refchannel = 'all';
    opmeeg_cleaned = ft_preprocessing(cfg,opmeeg_cleaned);
elseif isfield(params,'eeg_reref')
    cfg = [];
    cfg.channel = 'all';
    cfg.reref = 'yes';
    cfg.method = 'bipolar';
    cfg.refchannel = {params.eeg_reref};
    opmeeg_cleaned = ft_preprocessing(cfg,opmeeg_cleaned);
end

% Append ECG/EOG channels
opmeeg_cleaned.label = vertcat(opmeeg_cleaned.label,ExG.label);
for i = 1:length(opmeeg_cleaned.trial)
    opmeeg_cleaned.trial{i} = vertcat(opmeeg_cleaned.trial{i}, ExG.trial{i}); 
end

smplinfo = opmeeg_cleaned.sampleinfo;

% Reject jump trials
cfg = [];
cfg.channel = {'EEG*'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg, badtrl_opmeeg_jump] = ft_badsegment(cfg, opmeeg_cleaned);
opmeeg_cleaned = ft_rejectartifact(cfg,opmeeg_cleaned);

badtrl_opmeeg_std = zeros(0,2);
% % Reject noisy trials
% cfg = [];
% cfg.channel = {'EEG*'};
% cfg.metric = 'std';
% cfg.threshold = params.eeg_std_threshold;
% [cfg, badtrl_opmeeg_std] = ft_badsegment(cfg, opmeeg_cleaned);
% opmeeg_cleaned = ft_rejectartifact(cfg,opmeeg_cleaned);

% Save bad channels 
save(fullfile(save_path, [params.sub '_opmeeg_badchs']), ...
    'badchs_opmeeg_flat', ...
    'badchs_opmeeg_neighbors', "-v7.3"); 

% Save bad trials
[~,idx]=ismember(opmeeg_cleaned.sampleinfo,badtrl_opmeeg_jump,'rows');
badtrl_opmeeg_jump = find(idx);
[~,idx]=ismember(opmeeg_cleaned.sampleinfo,badtrl_opmeeg_std,'rows');
badtrl_opmeeg_std = find(idx);
save(fullfile(save_path, [params.sub '_opmeeg_badtrls']), ...
    'badtrl_opmeeg_jump', ...
    'badtrl_opmeeg_std', "-v7.3"); 

%save(fullfile(save_path, [params.sub '_opm_cleaned']), 'opm_cleaned',"-v7.3");
%save(fullfile(save_path, [params.sub '_opmeeg_cleaned']), 'opmeeg_cleaned',"-v7.3"); disp('done');

%% Downsample
if isfield(params,'ds_freq') && ~isempty(params.ds_freq) && params.ds_freq~=1000
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    opm_cleaned = ft_resampledata(cfg, opm_cleaned);
    cfg = [];
    cfg.resamplefs = params.ds_freq;
    opmeeg_cleaned = ft_resampledata(cfg, opmeeg_cleaned);
end


end