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

n_triggers = length(opm_timelocked);

numdipoles = params.numdipoles;

%% Fit dipoles
for i_trigger = 1:n_triggers

    % MEG
    cfg = [];
    cfg.gridsearch      = 'yes';                      
    cfg.sourcemodel     = sourcemodel;           
    cfg.headmodel       = headmodel;    
    cfg.senstype        = 'meg';            
    cfg.channel         = 'megmag';         
    cfg.nonlinear       = 'yes';     
    cfg.numdipoles      = numdipoles;
    if numdipoles == 2
        cfg.symmetry    = 'x';
    end
    cfg.latency         = peak_squid{i_trigger}.peak_latency + [-0.01 0.01];
    cfg.dipfit.checkinside = 'yes';
    dipole = ft_dipolefitting(cfg, squid_timelocked{i_trigger});
    if numdipoles == 2
        if dipole.dip.pos(1,1) > 0 % first dipole is on the right -> flip order
            dipole.dip.pos([2 1],:) = dipole.dip.pos([1 2],:);
            dipole.dip.mom([4 5 6 1 2 3]) = dipole.dip.mom([1 2 3 4 5 6]);
        end
    end
    squidmag_dipole{i_trigger} = dipole;

    cfg.channel         = 'meggrad';           
    dipole = ft_dipolefitting(cfg, squid_timelocked{i_trigger});
    if numdipoles == 2
        if dipole.dip.pos(1,1) > 0 % first dipole is on the right -> flip order
            dipole.dip.pos([2 1],:) = dipole.dip.pos([1 2],:);
            dipole.dip.mom([4 5 6 1 2 3]) = dipole.dip.mom([1 2 3 4 5 6]);
        end
    end
    squidgrad_dipole{i_trigger} = dipole;

    % OPM
    cfg = [];
    cfg.gridsearch      = 'yes';                         
    cfg.sourcemodel     = sourcemodel;            
    cfg.headmodel       = headmodel;    
    cfg.senstype        = 'meg';            
    cfg.channel         = '*bz';        
    cfg.nonlinear       = 'yes';    
    cfg.numdipoles      = numdipoles;
    if numdipoles == 2
        cfg.symmetry    = 'x';
    end
    cfg.latency         = peak_opm{i_trigger}.peak_latency + [-0.01 0.01];   
    cfg.dipfit.checkinside = 'yes';
    dipole = ft_dipolefitting(cfg, opm_timelocked{i_trigger});
    if numdipoles == 2
        if dipole.dip.pos(1,1) > 0 % first dipole is on the right -> flip order
            dipole.dip.pos([2 1],:) = dipole.dip.pos([1 2],:);
            dipole.dip.mom([4 5 6 1 2 3]) = dipole.dip.mom([1 2 3 4 5 6]);
        end
    end
    opm_dipole{i_trigger} = dipole;

    % Plot OPM vs SQUID
    for i_dip = 1:numdipoles
        pos_mag(i_dip,:) = squidmag_dipole{i_trigger}.dip.pos(i_dip,:);
        [~,idx] = max(vecnorm(squidmag_dipole{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),:),2,1));
        ori_mag(i_dip,:) = squidmag_dipole{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),idx);

        pos_grad(i_dip,:) = squidgrad_dipole{i_trigger}.dip.pos(i_dip,:);
        [~,idx] = max(vecnorm(squidgrad_dipole{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),:),2,1));
        ori_grad(i_dip,:) = squidgrad_dipole{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),idx);

        pos_opm(i_dip,:) =opm_dipole{i_trigger}.dip.pos(i_dip,:);
        [~,idx] = max(vecnorm(opm_dipole{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),:),2,1));
        ori_opm(i_dip,:) = opm_dipole{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),idx);
    end

    if params.numdipoles == 1
        h = figure;
        ft_plot_dipole(pos_mag,ori_mag,'color',colors(1,:))
        hold on;
        ft_plot_dipole(pos_opm,ori_opm,'color',colors(2,:))
        ft_plot_dipole(pos_grad,ori_grad,'color',colors(3,:))
        ft_plot_headmodel(headmodel,'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm') 
        hold off
        title([params.trigger_labels{i_trigger} ' (SQMAG-OPM = ' num2str(norm(pos_mag-pos_opm)*10,'%.1f') 'mm / SQGRAD-OPM = ' num2str(norm(pos_grad-pos_opm)*10,'%.1f') 'mm)'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' peak_squid{i_trigger}.label '_dipfit_SQUIDvOPM_trig-' params.trigger_labels{i_trigger} '.jpg']))
        close
    elseif params.numdipoles == 2
        h = figure;
        ft_plot_dipole(pos_mag(1,:),ori_mag(1,:),'color',colors(1,:))
        hold on;
        ft_plot_dipole(pos_mag(2,:),ori_mag(2,:),'color',colors(1,:))
        ft_plot_dipole(pos_opm(1,:),ori_opm(1,:),'color',colors(2,:))
        ft_plot_dipole(pos_opm(2,:),ori_opm(2,:),'color',colors(2,:))
        ft_plot_dipole(pos_grad(1,:),ori_grad(1,:),'color',colors(3,:))
        ft_plot_dipole(pos_grad(2,:),ori_grad(2,:),'color',colors(3,:))
        ft_plot_headmodel(headmodel,'EdgeAlpha',0,'FaceAlpha',0.3,'FaceColor',[229 194 152]/256,'unit','cm') 
        hold off
        title([params.trigger_labels{i_trigger} ' (SQMAG-OPM = ' num2str(norm(pos_mag(1,:)-pos_opm(1,:))*10,'%.1f') 'mm / SQGRAD-OPM = ' num2str(norm(pos_grad(1,:)-pos_opm(1,:))*10,'%.1f') 'mm)'])
        saveas(h, fullfile(save_path, 'figs', [params.sub '_' peak_squid{i_trigger}.label '_dipfit_SQUIDvOPM_trig-' params.trigger_labels{i_trigger} '.jpg']))
        close
    end
end
close all

%% Plot phalanges jointly
% SQUID
params.modality = 'squidmag';
dip = squidmag_dipole;
pos = zeros(n_triggers,3);
ori= zeros(n_triggers,3);
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        pos(((i_trigger-1)*numdipoles)+i_dip,:) = dip{i_trigger}.dip.pos(i_dip,:);
        [~,idx] = max(vecnorm(dip{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),:),2,1));
        ori(((i_trigger-1)*numdipoles)+i_dip,:) = dip{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),idx);
    end
end
mean_pos = mean(pos,1);

h=figure;
h.Position = [100 100 1000 1000];
set(gca,'LooseInset',get(gca,'TightInset'));
%010
subplot(2,2,1)
hold on
tmp = pos;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'tag', 'y')
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
view(0,0)
axis equal; axis tight; axis vis3d; axis off
%100
subplot(2,2,2)
hold on
tmp = pos;
tmp(:,1) = mean_pos(1)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'tag', 'x')
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
view(90,0)
axis equal; axis tight; axis vis3d; axis off
%001
subplot(2,2,3)
hold on
tmp = pos;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'tag', 'z')       
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
axis equal; axis tight; axis vis3d; axis off
% text
subplot(2,2,4);
axis off; % Turn off axis
hold on
text(0, 0.88, 'SQUID-MAG', 'FontWeight', 'bold');
if numdipoles == 2
    text(0.28, 0.8, 'dip_L (mm)', 'FontWeight', 'bold');
    text(0.28 + 0.4, 0.8, 'dip_R (mm)', 'FontWeight', 'bold');
else
    text(0.28, 0.8, ['dip' num2str(i_dip) ' (mm)'], 'FontWeight', 'bold');
end
for i_trigger = 1:n_triggers
        text(0, 0.8-(i_trigger*0.05), [params.trigger_labels{i_trigger} ': '], 'FontWeight', 'bold'); 
    for i_dip = 1:numdipoles
        text(0.28 + (i_dip-1)*0.4, 0.8-(i_trigger*0.05), [num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,1),'%.1f') ' / ' num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,2),'%.1f') ' / ' num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,3),'%.1f')]); 
    end
end
hold off
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' peak_squid{i_trigger}.label '_dipfit_mri.jpg']))
close all

% OPM
params.modality = 'opm';
dip = opm_dipole;
pos = zeros(n_triggers,3);
ori= zeros(n_triggers,3);
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        pos(((i_trigger-1)*numdipoles)+i_dip,:) = dip{i_trigger}.dip.pos(i_dip,:);
        [~,idx] = max(vecnorm(dip{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),:),2,1));
        ori(((i_trigger-1)*numdipoles)+i_dip,:) = dip{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),idx);
    end
end
mean_pos = mean(pos,1);

h=figure;
h.Position = [100 100 1000 1000];
set(gca,'LooseInset',get(gca,'TightInset'));
%010
subplot(2,2,1)
hold on
tmp = pos;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'tag', 'y')
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
view(0,0)
axis equal; axis tight; axis vis3d; axis off
%100
subplot(2,2,2)
hold on
tmp = pos;
tmp(:,1) = mean_pos(1)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'tag', 'x')
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
view(90,0)
axis equal; axis tight; axis vis3d; axis off
%001
subplot(2,2,3)
hold on
tmp = pos;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'tag', 'z')       
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
axis equal; axis tight; axis vis3d; axis off
% text
subplot(2,2,4);
axis off; % Turn off axis
hold on
text(0, 0.88, 'OPM', 'FontWeight', 'bold');
if numdipoles == 2
    text(0.28, 0.8, 'dip_L (mm)', 'FontWeight', 'bold');
    text(0.28 + 0.4, 0.8, 'dip_R (mm)', 'FontWeight', 'bold');
else
    text(0.28, 0.8, ['dip' num2str(i_dip) ' (mm)'], 'FontWeight', 'bold');
end
for i_trigger = 1:n_triggers
        text(0, 0.8-(i_trigger*0.05), [params.trigger_labels{i_trigger} ': '], 'FontWeight', 'bold'); 
    for i_dip = 1:numdipoles
        text(0.28 + (i_dip-1)*0.4, 0.8-(i_trigger*0.05), [num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,1),'%.1f') ' / ' num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,2),'%.1f') ' / ' num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,3),'%.1f')]); 
    end
end
hold off
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' peak_opm{i_trigger}.label '_dipfit_mri.jpg']))
close all

% SQUID-GRAD
params.modality = 'squidgrad';
dip = squidgrad_dipole;
pos = zeros(n_triggers,3);
ori= zeros(n_triggers,3);
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        pos(((i_trigger-1)*numdipoles)+i_dip,:) = dip{i_trigger}.dip.pos(i_dip,:);
        [~,idx] = max(vecnorm(dip{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),:),2,1));
        ori(((i_trigger-1)*numdipoles)+i_dip,:) = dip{i_trigger}.dip.mom((1:3)+((i_dip-1)*3),idx);
    end
end
mean_pos = mean(pos,1);

h=figure;
h.Position = [100 100 1000 1000];
set(gca,'LooseInset',get(gca,'TightInset'));
%010
subplot(2,2,1)
hold on
tmp = pos;
tmp(:,2) = mean_pos(2)-1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 1 0], 'tag', 'y')
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
view(0,0)
axis equal; axis tight; axis vis3d; axis off
%100
subplot(2,2,2)
hold on
tmp = pos;
tmp(:,1) = mean_pos(1)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [1 0 0], 'tag', 'x')
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
view(90,0)
axis equal; axis tight; axis vis3d; axis off
%001
subplot(2,2,3)
hold on
tmp = pos;
tmp(:,3) = mean_pos(3)+1;
ft_plot_slice(mri.anatomy, 'transform', mri.transform, 'location', mean_pos, 'orientation', [0 0 1], 'tag', 'z')       
for i_trigger = 1:n_triggers
    for i_dip = 1:numdipoles
        ft_plot_dipole(tmp((numdipoles*(i_trigger-1))+i_dip,:),ori((numdipoles*(i_trigger-1))+i_dip,:),'color',colors(i_trigger,:))
    end
end
hold off
axis equal; axis tight; axis vis3d; axis off
% text
subplot(2,2,4);
axis off; % Turn off axis
hold on
text(0, 0.88, 'SQUID-GRAD', 'FontWeight', 'bold');
if numdipoles == 2
    text(0.28, 0.8, 'dip_L (mm)', 'FontWeight', 'bold');
    text(0.28 + 0.4, 0.8, 'dip_R (mm)', 'FontWeight', 'bold');
else
    text(0.28, 0.8, ['dip' num2str(i_dip) ' (mm)'], 'FontWeight', 'bold');
end
for i_trigger = 1:n_triggers
        text(0, 0.8-(i_trigger*0.05), [params.trigger_labels{i_trigger} ': '], 'FontWeight', 'bold'); 
    for i_dip = 1:numdipoles
        text(0.28 + (i_dip-1)*0.4, 0.8-(i_trigger*0.05), [num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,1),'%.1f') ' / ' num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,2),'%.1f') ' / ' num2str(10*pos(((i_trigger-1)*numdipoles)+i_dip,3),'%.1f')]); 
    end
end
hold off
saveas(h, fullfile(save_path, 'figs', [params.sub '_' params.modality '_' peak_squid{i_trigger}.label '_dipfit_mri.jpg']))
close all

%% Save 
dipoles = squidmag_dipole;
save(fullfile(save_path, [ 'squidmag_' peak_opm{i_trigger}.label '_dipoles']), 'dipoles'); disp('done');
dipoles = squidgrad_dipole;
save(fullfile(save_path, [ 'squidgrad_' peak_opm{i_trigger}.label '_dipoles']), 'dipoles'); disp('done');
dipoles = opm_dipole;
save(fullfile(save_path, [ 'opm_' peak_opm{i_trigger}.label '_dipoles']),'dipoles'); disp('done');


end