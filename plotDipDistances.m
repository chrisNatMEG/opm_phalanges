function plotDipDistances(data,label,params,save_path,showSig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dat = median(data,3,'omitnan');

h = figure;
imagesc(dat);
colormap('winter');
colorbar;
axis square;

n_triggers = length(params.trigger_labels);

for i = 1:n_triggers
    for j = 1:n_triggers
        [~, p] = ttest(squeeze(data(i,j,:))); 
        txt = sprintf('%.2f',dat(i,j));
        if p<0.05 && showSig
            txt = [txt '*'];
        end
        text(j,i, txt,...
            'HorizontalAlignment','center','Color','w');
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