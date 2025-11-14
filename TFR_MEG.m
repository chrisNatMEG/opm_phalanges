function TFR_MEG(timelocked, save_path, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
n_triggers = length(params.trigger_codes);
TFR = cell(n_triggers,1);

for i_trigger = 1:n_triggers
    cfg = [];
    cfg.channel    = params.chs;
    cfg.method     = 'mtmfft';
    cfg.taper      = 'hanning';
    cfg.pad = 2;
    cfg.foi        = 34:1:46;         
    TFR{i_trigger} = ft_freqanalysis(cfg, timelocked{i_trigger});    % visual stimuli
    TFR{i_trigger}.SNR = max(max(TFR{i_trigger}.powspctrm)) / (mean(mean(TFR{i_trigger}.powspctrm(:,[1 end]))));
    TFR{i_trigger}.maxpow = max(max(TFR{i_trigger}.powspctrm));

    h=figure;
    plot(TFR{i_trigger}.freq,TFR{i_trigger}.powspctrm);
    xlabel('Freq [Hz]')
    ylabel('Power [T^2/Hz]')
    title(['FFT: ' params.trigger_labels{i_trig} ' (SNR=' num2str(max(max(TFR{i_trigger}.powspctrm)) /(mean(mean(TFR{i_trigger}.powspctrm(:,[1 end])))),'%.1f') ')' ])
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_FFT_trig-' params.trigger_labels{i_trig} '.jpg']))
    close all

    cfg = [];
    cfg.layout       = params.layout;
    h=figure; ft_topoplotER(cfg, TFR{i_trigger});
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_FFT-topo_trig-' params.trigger_labels{i_trig} '.jpg']))
    close all
end

% Save data
save(fullfile(save_path, [params.sub '_' params.modality '_TFR']), 'TFR', '-v7.3'); 
end