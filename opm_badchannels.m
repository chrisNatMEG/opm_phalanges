function [badchs, badchs_flat, badchs_std, badchs_neighbors, badchs_outlier, badchs_noloc] = opm_badchannels(data, trl, params, save_path)
%opm_badchannels Detects channels that are flat, have low correlation with
%thei_chs_gradeighbors or show a lot of jumping artifacts.
%   cfg.z_threshold
%   cfg.corr_threshold
%   cfg.n_neighbors
%   cfg.njump_threshold

std_threshold = ft_getopt(params, 'std_threshold', 5e-9);
n_neighbors     = ft_getopt(params, 'n_neighbors', 4);
corr_threshold  = ft_getopt(params, 'corr_threshold', 0.6);

cfg = [];
cfg.channel = '*_b*';
data = ft_selectdata(cfg, data);

chs = find(contains(data.label,'_b'));

%% Find channels with flat segments or high std
cfg = [];
cfg.length = 1;
data_seg = ft_redefinetrial(cfg,data);
for i_trl = 1:length(data_seg.trial)
    trl_std(:,i_trl) = std(data_seg.trial{i_trl},0,2);
end
badchs_flat = find(any(trl_std<1e-15,2));
badchs_std = find(mean(trl_std,2)>std_threshold);

badchs_noloc = find(startsWith(data.label(chs),'s'));

%% Neighbors
goodchs = setdiff(chs,[badchs_flat; badchs_std; badchs_noloc]);

cfg = [];
cfg.resamplefs = 200;
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;
cfg.lpinstabilityfix  = 'reduce';
data_lp = ft_resampledata(cfg,data);

cfg = [];
cfg.length = 1;
data_lp = ft_redefinetrial(cfg,data_lp);

% Create neighbor structure
goodchsZ = goodchs(contains(data_lp.label(goodchs),'bz'));
n_chs = length(goodchsZ);
[~,~,ib] = intersect(data_lp.label(goodchsZ),data_lp.grad.label,'stable');
chanpos = data_lp.grad.chanpos(ib,:);
neighbors = zeros(n_chs,n_neighbors);
for i = 1:size(chanpos,1)
        [~,tmp]= sort(vecnorm(chanpos-repmat(chanpos(i,:),[length(goodchsZ) 1]),2,2));
        neighbors(i,:) = tmp(2:(n_neighbors+1));
end
neighborscorr = zeros(n_chs,n_neighbors,length(data_lp.trial));
for i_trl = 1:length(data_lp.trial)
    dat = data_lp.trial{i_trl}(goodchsZ,:);
    for i = 1:n_chs
        for j = 1:n_neighbors
            tmp2 = corrcoef(dat(i,:),dat(int32(neighbors(i,j)),:));
            neighborscorr(i,j,i_trl) = abs(tmp2(1,2));
        end
    end 
end
badchs_neighbors = goodchsZ(max(mean(neighborscorr,3),[],2)<corr_threshold); % bad if no neighbors exceed correlation threshold

tmp = [];
for i = 1:length(badchs_neighbors)
    tmp = [tmp; find(contains(data_lp.label,data_lp.label{badchs_neighbors(i)}(1:end-2)))];
end
badchs_neighbors = tmp;
badchs = [badchs_flat; badchs_std; badchs_noloc; badchs_neighbors]; 

%% Epoch
cfg = [];
cfg.trl = trl;
data_epo = ft_redefinetrial(cfg,data);

cfg = [];
cfg.demean = 'yes';
cfg.dftfilter       = 'yes';        
cfg.dftfreq         = [40 50 60 80 100 ];
data_epo = ft_preprocessing(cfg,data_epo);

cfg = [];
cfg.channel = '*_b*';
cfg.metric = 'maxzvalue';
cfg.preproc.medianfilter  = 'yes';
cfg.preproc.medianfiltord  = 9;
cfg.preproc.absdiff       = 'yes';
cfg.threshold = params.z_threshold;
[cfg, ~] = ft_badsegment(cfg, data_epo);
data_epo = ft_rejectartifact(cfg,data_epo);

%% Spectrum
goodchs = setdiff(chs,badchs);
goodchsX = goodchs(contains(data_epo.label(goodchs),'bx'));
goodchsY = goodchs(contains(data_epo.label(goodchs),'by'));
goodchsZ = goodchs(contains(data_epo.label(goodchs),'bz'));
badchsX = badchs(contains(data_epo.label(badchs),'bx'));
badchsY = badchs(contains(data_epo.label(badchs),'by'));
badchsZ = badchs(contains(data_epo.label(badchs),'bz'));

cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foilim = [3 70];
freq = ft_freqanalysis(cfg, data_epo);
freq.powspctrm = sqrt(freq.powspctrm);

if ~isempty(goodchsX)
    badchs_outlierX = [];
    for it = 1:5
        tmp = freq.powspctrm(goodchsX,:);
        zscore = abs((tmp-mean(tmp,1))./std(tmp,0,1));
        badchs_outlierX = [badchs_outlierX; goodchsX(mean(zscore>params.outlier_zscore,2)>params.outlier_ratio)];
        goodchsX = setdiff(goodchsX,badchs_outlierX);
        if ~any(mean(zscore>params.outlier_zscore,2)>params.outlier_ratio)
            disp(['Outlier converged on iteration ' num2str(it)])
            break;
        elseif (it==5)
            disp('Outlier did not converge!!')
        end
    end
    thresholdX = mean(tmp,1) + params.outlier_zscore*std(tmp,0,1);
    h = figure;
    semilogy(freq.freq,freq.powspctrm(goodchsX,:),'Color',[211 211 211]/255)
    hold on
    if ~isempty(badchs_outlierX)
        semilogy(freq.freq,freq.powspctrm(badchs_outlierX,:),'r')
    end
    if ~isempty(badchsX)
        semilogy(freq.freq,freq.powspctrm(badchsX,:),'Color',[1 0.5 0])
    end
    semilogy(freq.freq,thresholdX,'k:')
    hold off
    title('X')
    saveas(h, fullfile(save_path, 'figs', [params.paradigm '_badchs_outliers_bX.jpg']))
    close all
else
    badchs_outlierX = [];
end

if ~isempty(goodchsY)
    badchs_outlierY = [];
    for it = 1:5
        tmp = freq.powspctrm(goodchsY,:);
        zscore = abs((tmp-mean(tmp,1))./std(tmp,0,1));
        badchs_outlierY = [badchs_outlierY; goodchsY(mean(zscore>params.outlier_zscore,2)>params.outlier_ratio)];
        goodchsY = setdiff(goodchsY,badchs_outlierY);
        if ~any(mean(zscore>params.outlier_zscore,2)>params.outlier_ratio)
            disp(['Outlier converged on iteration ' num2str(it)])
            break;
        elseif (it==5)
            disp('Outlier did not converge!!')
        end
    end
    thresholdY = mean(tmp,1) + params.outlier_zscore*std(tmp,0,1);
    h = figure;
    semilogy(freq.freq,freq.powspctrm(goodchsY,:),'Color',[211 211 211]/255)
    hold on
    if ~isempty(badchs_outlierY)
        semilogy(freq.freq,freq.powspctrm(badchs_outlierY,:),'r')
    end
    if ~isempty(badchsY)
        semilogy(freq.freq,freq.powspctrm(badchsY,:),'Color',[1 0.5 0])
    end
    semilogy(freq.freq,thresholdY,'k:')
    hold off
    title('Y')
    saveas(h, fullfile(save_path, 'figs', [params.paradigm '_badchs_outliers_bY.jpg']))
    close all
else
    badchs_outlierY = [];
end

if ~isempty(goodchsZ)
    badchs_outlierZ = [];
    for it = 1:5
        tmp = freq.powspctrm(goodchsZ,:);
        zscore = abs((tmp-mean(tmp,1))./std(tmp,0,1));
        badchs_outlierZ = [badchs_outlierZ; goodchsZ(mean(zscore>params.outlier_zscore,2)>params.outlier_ratio)];
        goodchsZ = setdiff(goodchsZ,badchs_outlierZ);
        if ~any(mean(zscore>params.outlier_zscore,2)>params.outlier_ratio)
            disp(['Outlier converged on iteration ' num2str(it)])
            break;
        elseif (it==5)
            disp('Outlier did not converge!!')
        end
    end
    thresholdZ = mean(tmp,1) + params.outlier_zscore*std(tmp,0,1);
    h = figure;
    semilogy(freq.freq,freq.powspctrm(goodchsZ,:),'Color',[211 211 211]/255)
    hold on
    if ~isempty(badchs_outlierZ)
        semilogy(freq.freq,freq.powspctrm(badchs_outlierZ,:),'r')
    end
    if ~isempty(badchsZ)
        semilogy(freq.freq,freq.powspctrm(badchsZ,:),'Color',[1 0.5 0])
    end
    semilogy(freq.freq,thresholdZ,'k:')
    hold off
    title('Z')
    saveas(h, fullfile(save_path, 'figs', [params.sub  '_opm_badchs_outliers_bZ.jpg']))
    close all
else
    badchs_outlierZ = [];
end
badchs_outlier = [badchs_outlierX; badchs_outlierY; badchs_outlierZ];

badchs = [badchs_flat; badchs_std; badchs_outlier; badchs_noloc; badchs_neighbors]; 

% Convert to channel labels
badchs = data.label(chs(badchs));
badchs_flat = data.label(chs(badchs_flat));
badchs_std = data.label(chs(badchs_std));
badchs_neighbors = data.label(chs(badchs_neighbors));
badchs_noloc = data.label(chs(badchs_noloc));
badchs_outlier = data.label(chs(badchs_outlier));