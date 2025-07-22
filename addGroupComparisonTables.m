function rptTable = addGroupComparisonTables(dataCells, triggerLabels, modalityLabels, formatStr)
    if nargin < 4
        formatStr = '%.2f';
    end

    import mlreportgen.report.*
    import mlreportgen.dom.*

    numModalities = numel(dataCells);
    numTriggers = numel(triggerLabels);
    numParticipants = size(dataCells{1}, 1);

    % Generate tables for each modality
    for i = 1:numModalities
        data = dataCells{i};
        meanOverParticipants = mean(data, 1);
        meanOverTriggers = mean(data, 2);
        grandMean = mean(data(:));

        % Create table with extra row and column for means
        tableData = cell(numParticipants + 1, numTriggers + 1);
        for r = 1:numParticipants
            for c = 1:numTriggers
                tableData{r, c} = sprintf(formatStr, data(r, c));
            end
            tableData{r, end} = sprintf(formatStr, meanOverTriggers(r));
        end
        for c = 1:numTriggers
            tableData{end, c} = sprintf(formatStr, meanOverParticipants(c));
        end
        tableData{end, end} = sprintf(formatStr, grandMean);

        % Add header row and column
        headerRow = [modalityLabels{i}, triggerLabels, {'Mean'}];
        tableData = [headerRow; [strcat("sub_", string(1:numParticipants))'; "Mean"], tableData];

        rptTable{i} = FormalTable(tableData);
        %rptTable.Title = sprintf('Modality: %s', modalityLabels{i});
        rptTable{i}.Header.Style = {Bold(true)};
        rptTable{i}.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
    end

    % Create comparison table with significance stars
    comparisonData = zeros(numModalities, numTriggers+1);
    for i = 1:numModalities
        comparisonData(i, 1:(end-1)) = mean(dataCells{i}, 1);
        comparisonData(i, end) = mean(mean(dataCells{i}, 1));
    end

    comparisonTable = cell(numModalities + 1, numTriggers + 1 +1);
    comparisonTable(1, 2:end-1) = triggerLabels;
    comparisonTable(1, end) = {'Mean'};
    comparisonTable(2:end, 1) = modalityLabels(:);

    for r = 1:numModalities
        for c = 1:(numTriggers+1)
            val = comparisonData(r, c);
            comparisonTable{r+1, c+1} = sprintf(formatStr, val);
        end
    end

    % Add significance stars
    for c = 1:numTriggers
        for i = 1:numModalities
            stars = "";
            for j = 1:numModalities
                if i ~= j
                    [~, p] = ttest(dataCells{i}(:, c), dataCells{j}(:, c));
                    if p < 0.001
                        stars = stars + "***";
                    elseif p < 0.01
                        stars = stars + "**";
                    elseif p < 0.05
                        stars = stars + "*";
                    else
                        stars = stars + " ";
                    end
                    stars = stars + "/";
                end
            end
            if ~isempty(stars)
                comparisonTable{i+1, c+1} = sprintf('%s %s', comparisonTable{i+1, c+1}, stars{1}(1:end-1));
            end
        end
    end
    for i = 1:numModalities
        stars = "";
        for j = 1:numModalities
            if i ~= j
                [~, p] = ttest(mean(dataCells{i}(:, :),2), mean(dataCells{j}(:, :),2));
                if p < 0.001
                    stars = stars + "***";
                elseif p < 0.01
                    stars = stars + "**";
                elseif p < 0.05
                    stars = stars + "*";
                else
                    stars = stars + " ";
                end
                if j < numModalities
                    stars = stars + "/";
                end
            end
        end
        if ~isempty(stars)
            comparisonTable{i+1, c+2} = sprintf('%s %s', comparisonTable{i+1, c+2}, stars{1}(1:end-1));
        end
    end

    rptTable{i+1} = FormalTable(comparisonTable);
    %rptTable.Title = 'Trigger-wise Means with Significance';
    rptTable{i+1}.Header.Style = {Bold(true)};
    rptTable{i+1}.Style = {Border('double'), Width('100%'), RowSep('solid'), ColSep('solid')};
end
