function [squidmag_dipole, squidgrad_dipole, opm_dipole] = fit_dipoles(save_path,squid_timelocked,opm_timelocked,headmodel,mri,peak_squid, peak_opm, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

colors = [[0 0.4470 0.7410]; % blue
    [0.8500 0.3250 0.0980]; % red
    [0.9290 0.6940 0.1250]; % yellow
    [0.4940 0.1840 0.5560]; % purple
    [0.4660 0.6740 0.1880]; % green
    [0.6350 0.0780 0.1840]]; % light blue

cfg              = [];
cfg.resolution   = 0.5;
cfg.tight        = 'yes';
cfg.inwardshift  = 0;
cfg.headmodel    = headmodel;
sourcemodel    = ft_prepare_sourcemodel(cfg);

%% Fit dipoles
for i_phalange = 1:5
    opm_timelocked{i_phalange}.avg = -opm_timelocked{i_phalange}.avg;

    % MEG
    cfg = [];
    cfg.gridsearch      = 'yes';           
    cfg.numdipoles      = 1;                
    cfg.sourcemodel     = sourcemodel;           
    cfg.headmodel       = headmodel;    
    cfg.senstype        = 'meg';            
    cfg.channel         = 'megmag';         
    cfg.nonlinear       = 'yes';           
    cfg.latency         = peak_squid{i_phalange}.peak_latency + [-0.01 0.01];
    cfg.dipfit.checkinside = 'yes';
    squidmag_dipole{i_phalange} = ft_dipolefitting(cfg, squid_timelocked{i_phalange});
    
    cfg.latency         = peak_squid{i_phalange}.peak_latency + [-0.01 0.01];   
    cfg.channel         = 'meggrad';           
    squidgrad_dipole{i_phalange} = ft_dipolefitting(cfg, squid_timelocked{i_phalange});

    % OPM
    cfg = [];
    cfg.gridsearch      = 'yes';           
    cfg.numdipoles      = 1;                
    cfg.sourcemodel     = sourcemodel;            
    cfg.headmodel       = headmodel;    
    cfg.senstype        = 'meg';            
    cfg.channel         = '*bz';        
    cfg.nonlinear       = 'yes';            
    cfg.latency         = peak_opm{i_phalange}.peak_latency + [-0.01 0.01];   
    cfg.dipfit.checkinside = 'yes';
    opm_dipole{i_phalange} = ft_dipolefitting(cfg, opm_timelocked{i_phalange});

    % Plot OPM vs SQUID
    pos_mag = squidmag_dipole{i_phalange}.dip.pos;
    [~,idx] = max(vecnorm(squidmag_dipole{i_phalange}.dip.mom,2,1));
    ori_mag = squidmag_dipole{i_phalange}.dip.mom(:,idx);

    pos_gad = squidgrad_dipole{i_phalange}.dip.pos;
    [~,idx] = max(vecnorm(squidgrad_dipole{i_phalange}.dip.mom,2,1));
    ori_grad = squidgrad_dipole{i_phalange}.dip.mom(:,idx);
    
    pos_opm = opm_dipole{i_phalange}.dip.pos;
    %pos_opm = opm_trans.transformPointsInverse(pos_opm);
    [~,idx] = max(vecnorm(opm_dipole{i_phalange}.dip.mom,2,1));
    ori_opm = -opm_dipole{i_phalange}.dip.mom(:,idx);

    h = figure;
    ft_plot_dipole(pos_mag,ori_mag,'color',colors(1,:))
    hold on;
    ft_plot_dipole(pos_opm,ori_opm,'color',colors(2,:))
    ft_plot_dipole(pos_gad,ori_grad,'color',colors(3,:))
    ft_plot_headmodel(headmodel,'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm') 
    hold off
    title([params.phalange_labels{i_phalange} ' (SQMAG-OPM = ' num2str(norm(pos_mag-pos_opm)*10,'%.1f') 'mm / SQGRAD-OPM = ' num2str(norm(pos_gad-pos_opm)*10,'%.1f') 'mm)'])
    %legend('SQUIDMAG','OPM','SQUIDPLANAR','brain')
    saveas(h, fullfile(save_path, 'figs', [params.sub '_' peak_squid{i_phalange}.label '_dipfit_SQUIDvOPM_ph-' params.phalange_labels{i_phalange} '.jpg']))
    close
end
close all

%% Plot phalanges jointly
% SQUID
params.modality = 'squidmag';
pos_mag = zeros(5,3);
ori_mag = zeros(5,3);
for i = 1:5
    pos_mag(i,:) = squidmag_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(squidmag_dipole{i}.dip.mom,2,1));
    ori_mag(i,:) = squidmag_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_mag,1);

h=figure;
h.Position = [100 100 800 800];
set(gca,'LooseInset',get(gca,'TightInset'));
%010
subplot(2,2,1)
hold on
tmp = pos_mag;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'tag', 'y')
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
view(0,0)
axis equal; axis tight; axis vis3d; axis off
%100
subplot(2,2,2)
hold on
tmp = pos_mag;
tmp(:,1) = mean_pos(1)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'tag', 'x')
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
view(90,0)
axis equal; axis tight; axis vis3d; axis off
%001
subplot(2,2,3)
hold on
tmp = pos_mag;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'tag', 'z')       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_mag(i,:),'color',colors(i,:))
end
hold off
axis equal; axis tight; axis vis3d; axis off
h.Position = [100 100 800 900];
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' peak_squid{i_phalange}.label '_dipfit_mri.jpg']))
close all

% OPM
params.modality = 'opm';
pos_opm = zeros(5,3);
ori_opm = zeros(5,3);
for i = 1:5
    pos_opm(i,:) = opm_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(opm_dipole{i}.dip.mom,2,1));
    ori_opm(i,:) = -opm_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_opm,1);

h=figure;
h.Position = [100 100 800 800];
set(gca,'LooseInset',get(gca,'TightInset'));
%010
subplot(2,2,1)
hold on
tmp = pos_opm;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'tag', 'y')
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
view(0,0)
axis equal; axis tight; axis vis3d; axis off
%100
subplot(2,2,2)
hold on
tmp = pos_opm;
tmp(:,1) = mean_pos(1)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'tag', 'x')
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
view(90,0)
axis equal; axis tight; axis vis3d; axis off
%001
subplot(2,2,3)
hold on
tmp = pos_opm;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'tag', 'z')       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_opm(i,:),'color',colors(i,:))
end
hold off
axis equal; axis tight; axis vis3d; axis off
h.Position = [100 100 800 900];
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' peak_opm{i_phalange}.label '_dipfit_mri.jpg']))
close all

% SQUID-GRAD
params.modality = 'squidgrad';
pos_grad = zeros(5,3);
ori_grad = zeros(5,3);
for i = 1:5
    pos_grad(i,:) = squidgrad_dipole{i}.dip.pos;
    [~,idx] = max(vecnorm(squidgrad_dipole{i}.dip.mom,2,1));
    ori_grad(i,:) = squidgrad_dipole{i}.dip.mom(:,idx);
end
mean_pos = mean(pos_grad,1);

h=figure;
h.Position = [100 100 800 800];
set(gca,'LooseInset',get(gca,'TightInset'));
%010
subplot(2,2,1)
hold on
tmp = pos_grad;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'tag', 'y')
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_grad(i,:),'color',colors(i,:))
end
hold off
view(0,0)
axis equal; axis tight; axis vis3d; axis off
%100
subplot(2,2,2)
hold on
tmp = pos_grad;
tmp(:,1) = mean_pos(1)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'tag', 'x')
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_grad(i,:),'color',colors(i,:))
end
hold off
view(90,0)
axis equal; axis tight; axis vis3d; axis off
%001
subplot(2,2,3)
hold on
tmp = pos_grad;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'tag', 'z')       
for i = 1:5
    ft_plot_dipole(tmp(i,:),ori_grad(i,:),'color',colors(i,:))
end
hold off
axis equal; axis tight; axis vis3d; axis off
h.Position = [100 100 800 900];
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' peak_squid{i_phalange}.label '_dipfit_mri.jpg']))
close all

%% Save 
save(fullfile(save_path, [ peak_squid{i_phalange}.label '_dipoles']), 'squidmag_dipole', 'squidgrad_dipole', 'opm_dipole'); disp('done');

end