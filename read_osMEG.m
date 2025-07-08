function [opm_cleaned, opmeeg_cleaned] = read_osMEG(opm_file, aux_file, save_path, params)
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
trl_opm(:,1:2) = trl_opm(:,1:2) + floor(0.041*opm_raw.fsample); % adjust for stim delay
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
trl_aux(:,1:2) = trl_aux(:,1:2) + floor(0.041*aux_raw.fsample); % adjust for stim delay
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
cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
if ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
    if params.filter.hp_freq<1
        cfg.hpfilttype = 'firws';
    end
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

cfg.lpfilter        = 'yes';         
cfg.lpfreq          = params.filter.lp_freq;
if ~isempty(params.filter.hp_freq)
    cfg.hpfilter        = 'yes'; 
    cfg.hpfreq          = params.filter.hp_freq;
    cfg.hpinstabilityfix  = 'reduce';
%     if params.filter.hp_freq<1
%         cfg.hpfilttype = 'firws';
%     end
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
cfg = [];
cfg.trl = trl_opm;
cfg.z_threshold = params.z_threshold;
cfg.corr_threshold = params.corr_threshold;
[badchs_opm, badchs_opm_flat, badchs_opm_std, badchs_opm_neighbors, badchs_opm_zmax, badchs_outlier, badtrl_opm_zmax] = opm_badchannels(cfg,opm_raw);
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

%% Flip? 
chs = find(contains(comb.label,'bz'));
if isfield(params,'flip_sign') && params.flip_sign
    for i = 1:length(comb.trial)
        comb.trial{i}(chs,:) = -comb.trial{i}(chs,:);
    end
end

%% OPM 
cfg = [];
cfg.channel = setdiff(comb.label,badchs_opm);
cfg.trials  = setdiff(1:length(comb.trial),badtrl_opm_zmax); % remove bad trials
comb = ft_selectdata(cfg, comb);

tmp = comb;
cfg = [];
cfg.demean          = 'yes';
cfg.baselinewindow  = [-params.pre 0];
tmp = ft_preprocessing(cfg,tmp);
cfg = [];
cfg.channel = '*bz';
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [1 100];
freq = ft_freqanalysis(cfg, tmp);
clear tmp
h = figure;
semilogy(freq.freq,freq.powspctrm)
xlabel('Frequency (Hz)')
ylabel('Power (T^2)')
title('OPM spectrum - pre')
saveas(h, fullfile(save_path, 'figs', [params.sub '_opm_spectrum_filt.jpg']))
close all

% HFC
% i_chs = find(contains(op_cleaned.label,'bz'));
% chs = opm_cleaned.label(i_chs);
% ori = zeros(length(chs),3);
% for i = 1:length(chs)
%     ori(i,:) = opm_cleaned.grad.chanori(find(strcmp(opm_cleaned.grad.label,chs{i})),:);
%     i_chs_grad(i) = find(strcmp(opm_cleaned.grad.label,chs{i}));
% end
% opm_cleaned.grad.M = eye(size(ori,1)) - ori*pinv(ori);
% for i = 1:length(opm_cleaned.trial)
%     opm_cleaned.trial{i}(i_chs,:) = params.sign * opm_cleaned.grad.M*opm_cleaned.trial{i}(i_chs,:);
% end
% opm_cleaned.grad.tra(i_chs_grad,i_chs_grad) = opm_cleaned.grad.M * opm_cleaned.grad.tra(i_chs_grad,i_chs_grad); % update grad
cfg = []; % separate ExG channels
cfg.channel = {'EOG*', 'ECG*'};
ExG = ft_selectdata(cfg,comb);

%% Spatiotemporal filtering
if params.do_hfc
    cfg = [];
    cfg.channel = '*bz';
    cfg.order = params.hfc_order;
    cfg.residualcheck = 'no';
    opm_cleaned = ft_denoise_hfc(cfg, comb);
    
    tmp = opm_cleaned;
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-params.pre 0];
    tmp = ft_preprocessing(cfg,tmp);
    cfg = [];
    cfg.channel = '*bz';
    cfg.output = 'pow';
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.foilim = [1 100];
    freq = ft_freqanalysis(cfg, tmp);
    clear tmp
    h = figure;
    semilogy(freq.freq,freq.powspctrm)
    xlabel('Frequency (Hz)')
    ylabel('Power (T^2)')
    title('OPM spectrum - postHFC')
    saveas(h, fullfile(save_path, 'figs', [params.sub '_opm_spectrum_HFC.jpg']))
    close all

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
    
    tmp = opm_cleaned;
    cfg = [];
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-params.pre 0];
    tmp = ft_preprocessing(cfg,tmp);
    cfg = [];
    cfg.channel = '*bz';
    cfg.output = 'pow';
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.foilim = [1 100];
    freq = ft_freqanalysis(cfg, tmp);
    clear tmp
    h = figure;
    semilogy(freq.freq,freq.powspctrm)
    xlabel('Frequency (Hz)')
    ylabel('Power (T^2)')
    title('OPM spectrum - postAMM')
    saveas(h, fullfile(save_path, 'figs', [params.sub '_opm_spectrum_AMM.jpg']))
    close all
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

% Reject jump trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg,badtrl_opm_jump] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'std';
cfg.threshold = params.opm_std_threshold;
[cfg,badtrl_opm_std] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

% Reject noisy trials
cfg = [];
cfg.channel = {'*bz'};
cfg.metric = 'range';
cfg.threshold = params.opm_range_threshold;
[cfg,badtrl_opm_range] = ft_badsegment(cfg, opm_cleaned);
opm_cleaned = ft_rejectartifact(cfg,opm_cleaned);

% Convert grad unit to cm to match TRIUX grad
opm_cleaned.grad = ft_convert_units(opm_cleaned.grad,'cm');

%% EEG
cfg = [];
cfg.channel = 'EEG';
opmeeg_cleaned = ft_selectdata(cfg, comb);

% Interpolate bad chs
cfg = [];
cfg.method = 'triangulation';
cfg.senstype = 'EEG';
neighbors = ft_prepare_neighbours(cfg,opmeeg_cleaned);
cfg = [];
cfg.method = 'spline';
cfg.neighbors = neighbors;
cfg.badchannel = badchs_opmeeg;
cfg.senstype = 'EEG';
opmeeg_cleaned = ft_channelrepair(cfg, opmeeg_cleaned);

% Re-reference
cfg = [];
cfg.channel = 'all';
cfg.refef = 'yes';
cfg.refmethod = 'avg';
cfg.reffchannel = 'all';
opmeeg_cleaned = ft_preprocessing(cfg,opmeeg_cleaned);

% Append ECG/EOG channels
opmeeg_cleaned.label = vertcat(opmeeg_cleaned.label,ExG.label);
for i = 1:length(opmeeg_cleaned.trial)
    opmeeg_cleaned.trial{i} = vertcat(opmeeg_cleaned.trial{i}, ExG.trial{i}); 
end

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

% Reject noisy trials
cfg = [];
cfg.channel = {'EEG*'};
cfg.metric = 'std';
cfg.threshold = params.eeg_std_threshold;
[cfg, badtrl_opmeeg_std] = ft_badsegment(cfg, opmeeg_cleaned);
opmeeg_cleaned = ft_rejectartifact(cfg,opmeeg_cleaned);

%% Spectra
cfg = [];
cfg.channel = 'EEG*';
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [1 100];
freq = ft_freqanalysis(cfg, opmeeg_cleaned);
h = figure;
semilogy(freq.freq,freq.powspctrm)
xlabel('Frequency (Hz)')
ylabel('Power (T^2)')
title('OPM-EEG spectrum - preICA')
saveas(h, fullfile(save_path, 'figs', [params.sub '_opmeeg_spectrum_preICA.jpg']))
close all

cfg = [];
cfg.channel = '*bz';
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [1 100];
freq = ft_freqanalysis(cfg, opm_cleaned);
h = figure;
semilogy(freq.freq,freq.powspctrm)
xlabel('Frequency (Hz)')
ylabel('Power (T^2)')
title('OPM spectrum - preICA')
saveas(h, fullfile(save_path, 'figs', [params.sub '_opm_spectrum_preICA.jpg']))
close all

%% Save 
save(fullfile(save_path, [params.sub '_opm_badchs']), ...
    'badchs_opm_flat', ...
    'badchs_opm_std', ...
    'badchs_opm_neighbors', ...
    'badchs_opm_zmax' , ...
    'badchs_outlier',"-v7.3"); 

save(fullfile(save_path, [params.sub '_opmeeg_badchs']), ...
    'badchs_opmeeg_flat', ...
    'badchs_opmeeg_neighbors', "-v7.3"); 

[~,idx]=ismember(opm_cleaned.sampleinfo,badtrl_opm_jump,'rows');
badtrl_opm_jump = find(idx);
[~,idx]=ismember(opm_cleaned.sampleinfo,badtrl_opm_std,'rows');
badtrl_opm_std = find(idx);
[~,idx]=ismember(opm_cleaned.sampleinfo,badtrl_opm_range,'rows');
badtrl_opm_range = find(idx);
save(fullfile(save_path, [params.sub '_opm_badtrls']), ...
    'badtrl_opm_jump', ...
    'badtrl_opm_std', ...
    'badtrl_opm_range', ...
    'badtrl_opm_zmax' ,"-v7.3"); 

[~,idx]=ismember(opmeeg_cleaned.sampleinfo,badtrl_opmeeg_jump,'rows');
badtrl_opmeeg_jump = find(idx);
[~,idx]=ismember(opmeeg_cleaned.sampleinfo,badtrl_opmeeg_std,'rows');
badtrl_opmeeg_std = find(idx);
save(fullfile(save_path, [params.sub '_opmeeg_badtrls']), ...
    'badtrl_opmeeg_jump', ...
    'badtrl_opmeeg_std', "-v7.3"); 

%save(fullfile(save_path, [params.sub '_opm_cleaned']), 'opm_cleaned',"-v7.3");
%save(fullfile(save_path, [params.sub '_opmeeg_cleaned']), 'opmeeg_cleaned',"-v7.3"); disp('done');

end