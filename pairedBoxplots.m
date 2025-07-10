function pairedBoxplots(data, triggerLabels, yLabelStr, titleStr, save_path, showSig)
    % data: cell array containing 2 or 3 matrices (participants x triggercodes)
    % triggerLabels: cell array of strings for x-axis labels
    % yLabelStr: string for y-axis label
    % titleStr: string for plot title

    % Validate input
    nGroups = numel(data);
    if nGroups < 2 || nGroups > 3
        error('Function supports comparison of 2 or 3 matrices.');
    end

    [nParticipants, nTriggers] = size(data{1});
    if length(triggerLabels) ~= nTriggers
        error('Number of trigger labels must match number of trigger codes.');
    end

    % Prepare data for boxplot
    allData = [];
    group = [];
    position = [];
    offset = linspace(-0.2, 0.2, nGroups);  % spacing between groups

    for g = 1:nGroups
        for t = 1:nTriggers
            idx = (t - 1) * nParticipants + (1:nParticipants);
            allData(idx, 1) = data{g}(:, t);
            group(idx, 1) = t + offset(g);
            position(idx, 1) = t;
        end
    end

    % Plot
    h = figure;
    hold on;
    colors = lines(nGroups);
    for g = 1:nGroups
        for t = 1:nTriggers
            x = t + offset(g);
            boxplot(data{g}(:, t), 'positions', x, 'colors', colors(g,:), ...
                    'widths', 0.15);
        end
    end

    % Formatting
    set(gca, 'XTick', 1:nTriggers, 'XTickLabel', triggerLabels);
    ylabel(yLabelStr);
    title(titleStr);
    %legend(arrayfun(@(i) sprintf('Matrix %d', i), 1:nGroups, 'UniformOutput', false));
    xtickangle(45);
    grid on;

    if showSig
        meanOffset = mean(offset);
        sigPairs = {};
        pValues = [];

        %Compare each pair of groups
        for g1 = 1:nGroups-1
            for g2 = g1+1:nGroups
                for t = 1:nTriggers
                    [~, p] = ttest(data{g1}(:,t),data{g2}(:,t));
                    sigPairs{end+1} = [t + offset(g1), t+ offset(g2)];
                    pValues(end+1) = p;
                end
            end
        end

        % Add significance markers
        sigstar(sigPairs,pValues)
    end

    hold off;
    saveas(h, save_path);
    close
end
