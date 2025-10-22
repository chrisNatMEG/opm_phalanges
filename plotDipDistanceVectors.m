function plotDipDistanceVectors(data,label,params,save_path,showSig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dat = median(data,3,'omitnan'); % median over subs

h = figure;
imagesc(mean(dat,4));
colormap('winter');
colorbar;
axis square;

n_triggers = length(params.trigger_labels);
y_offsets = [-0.2, 0, 0.2];

for i = 1:n_triggers
    for j = 1:n_triggers
        p = [];
        [~, p(1)] = ttest(squeeze(data(i,j,:,1))); 
        [~, p(2)] = ttest(squeeze(data(i,j,:,2)));
        [~, p(3)] = ttest(squeeze(data(i,j,:,3)));
        for k = 1:3
            txt = [sprintf('%.1f',dat(i,j,1,k))];
            if p(k)<0.05 && showSig
                txt = [txt '*'];
            end
            text(j,i + y_offsets(k), txt,...
                'HorizontalAlignment','center','Color','w');
        end
    end
end

set(gca, 'XTick', 1:n_triggers, 'XTickLabel', params.trigger_labels)
set(gca, 'YTick', 1:n_triggers, 'YTickLabel', params.trigger_labels)
xlabel('Trigger')
ylabel('Trigger')
title(['Median distance between phalanges - ' label])
saveas(h, save_path);

close all
end