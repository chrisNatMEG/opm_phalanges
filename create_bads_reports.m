function create_bads_reports(base_folder_path, subs, params)
% Import necessary packages
import mlreportgen.report.*
import mlreportgen.dom.*

% Loop over each subject in the list
for subNumber = subs
    % Convert subNumber to a two-digit string
    subStr = sprintf('%02d', subNumber);

    % Create the report name
    reportName = ['sub_', subStr, ' bads'];

    % Define the subject folder path
    subjectFolderPath = fullfile(base_folder_path, ['sub_', subStr]);

    % Create a new report
    rpt = Report(fullfile(subjectFolderPath, reportName), 'pdf');

    % Add a title page
    titlepg = TitlePage;
    titlepg.Title = ['Subject ', subStr, ' Preprocessing'];
    add(rpt, titlepg);

    % Add the table of contents
    toc = TableOfContents();
    add(rpt, toc);

    % Define sections
    sections = {'opm', 'meg', 'opmeeg', 'megeeg'};

    %% Bad Channels chapter
    chapter = Chapter('Bad Channels');
    chapter.Numbered = false; % Remove chapter numbering
    for i_section = 1:length(sections)
        filePath = fullfile(subjectFolderPath, ['sub_', subStr, '_' sections{i_section} '_badchs.mat']);
        if isfile(filePath)
            data = load(filePath);
            section = Section(sections(i_section));
            section.Numbered = false; % Remove section numbering
            fields = fieldnames(data);
            totalLength = 0;
            for j = 1:length(fields)
                varData = data.(fields{j});
                if isempty(varData)
                    varLength = 0;
                    varList = '';
                elseif iscell(varData)
                    varLength = length(varData);
                    varList = strjoin(varData, ', ');
                else
                    varLength = max(size(varData));
                    varList = num2str(varData(:)');
                end
                totalLength = totalLength + varLength;
                para = Paragraph([fields{j},' (n=', num2str(varLength), '): ', varList]);
                add(section, para);
            end
            para = Paragraph(['n_{total}: ', num2str(totalLength)]);
            add(section, para);
            add(chapter, section);
        end
    end
    add(rpt, chapter);

    %% Bad Trials chapter
    chapter = Chapter('Bad Trials');
    chapter.Numbered = false; % Remove chapter numbering
    for i_section = 1:length(sections)
        filePath = fullfile(subjectFolderPath, ['sub_', subStr, '_' sections{i_section} '_badtrls.mat']);
        if isfile(filePath)
            data = load(filePath);
            section = Section(sections(i_section));
            section.Numbered = false; % Remove section numbering
            fields = fieldnames(data);
            totalLength = 0;
            for j = 1:length(fields)
                varData = data.(fields{j});
                if isempty(varData)
                    varLength = 0;
                    varList = '';
                elseif iscell(varData)
                    varLength = length(varData);
                    varList = strjoin(varData, ', ');
                else
                    varLength = max(size(varData));
                    varList = strjoin(arrayfun(@(x) sprintf('%d %d', varData(x, 1), varData(x, 2)), 1:size(varData, 1), 'UniformOutput', false), '; ');
                end
                totalLength = totalLength + varLength;
                para = Paragraph([fields{j},' (n=', num2str(varLength), '): ', varList]);
                add(section, para);
            end
            para = Paragraph(['n_{total}: ', num2str(totalLength)]);
            add(section, para);
            add(chapter, section);
        end
    end
    add(rpt, chapter);

    %% ICA chapter
    chapter = Chapter('ICA');
    chapter.Numbered = false; % Remove chapter numbering
    for i_section = 1:length(sections)
        filePath = fullfile(subjectFolderPath,['sub_' subStr '_' sections{i_section} '_ica_comp.mat']);
        if isfile(filePath)
            data = load(filePath);
            section = Section(sections(i_section));
            section.Numbered = false; % Remove section numbering
            fields = {'ecg_comp_idx', 'eog1_comp_idx', 'eog2_comp_idx'};
            totalLength = 0;
            for j = 1:length(fields)
                varData = data.(fields{j});
                if isempty(varData)
                    varLength = 0;
                    varList = '';
                elseif iscell(varData)
                    varLength = length(varData);
                    varList = strjoin(varData, ', ');
                else
                    varLength = max(size(varData));
                    varList = num2str(varData(:)');
                    varList = strjoin(arrayfun(@num2str, varData(:)', 'UniformOutput', false), ', ');
                end
                totalLength = totalLength + varLength;
                para = Paragraph([fields{j}, ': ', varList]);
                add(section, para);
            end

            img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections{i_section} '_ica_rejected_comps1.jpg']));
            img.Style = {ScaleToFit};
            add(section, img);

            add(chapter, section);
        end
    end
    add(rpt, chapter);

    %% Timelocked chapter
    chapter = Chapter('Timelocked');
    chapter.Numbered = false; % Remove chapter numbering
    for i_section = 1:length(sections)
        section = Section(sections(i_section));
        section.Numbered = false; % Remove section numbering

        if contains(sections(i_section),'eeg')
            fieldMultiplier = 1e9;
            fieldUnit = '[nV]';
        else
            fieldMultiplier = 1e15;
            fieldUnit = '[fT]';
        end

        % Load the .mat file
        data = load(fullfile(subjectFolderPath,['sub_' subStr '_' sections{i_section} '_M100.mat']));
        M100 = data.M100;
        
        % Create the table with the required data
        num_phalanges = length(params.phalange_labels);
        T = table('Size', [6, num_phalanges], 'VariableTypes', repmat({'double'}, 1, num_phalanges), 'VariableNames', params.phalange_labels);
        T.Properties.RowNames = {['peak_amplitude ' fieldUnit], 'peak_latency [ms]', 'SNR_prestim', 'SNR_stderr', ['std_prestim ' fieldUnit], ['stderr ' fieldUnit]};
        
        for i_section = 1:num_phalanges
            phalange_data = M100{i_section};
            T{['peak_amplitude ' fieldUnit], params.phalange_labels{i_section}} = fieldMultiplier*phalange_data.max_amplitude;
            T{'peak_latency [ms]', params.phalange_labels{i_section}} = 1e3*phalange_data.peak_latency;
            T{'SNR_prestim', params.phalange_labels{i_section}} = phalange_data.max_amplitude / phalange_data.prestim_std;
            T{'SNR_stderr', params.phalange_labels{i_section}} = phalange_data.max_amplitude / phalange_data.std_error;
            T{['std_prestim ' fieldUnit], params.phalange_labels{i_section}} = fieldMultiplier*phalange_data.prestim_std;
            T{['stderr ' fieldUnit], params.phalange_labels{i_section}} = fieldMultiplier*phalange_data.std_error;
        end
        
        % Convert the MATLAB table to a DOM table
        domTable = Table();
        domTable.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
        
        % Add header row
        headerRow = TableRow();
        append(headerRow, TableEntry('Metric'));
        for i_section = 1:num_phalanges
            headerEntry = TableEntry(params.phalange_labels{i_section});
            headerEntry.Style = {HAlign('center'), Bold()};
            append(headerRow, headerEntry);
        end
        append(domTable, headerRow);
        
        % Add data rows
        for i_section = 1:height(T)
            row = TableRow();
            rowNameEntry = TableEntry(T.Properties.RowNames{i_section});
            rowNameEntry.Style = {Bold()};
            append(row, rowNameEntry);
            for j = 1:num_phalanges
                entry = TableEntry(num2str(T{i_section, j},'%.1f'));
                entry.Style = {HAlign('right')}; 
                append(row, entry);
            end
            append(domTable, row);
        end
        add(section,domTable);
        add(chapter,section)
    end
    add(rpt, chapter);

    %% Butterfly plots chapter
    chapter = Chapter('Butterfly plots');
    chapter.Numbered = false; % Remove chapter numbering
       
    for i_phalange = 1:length(params.phalange_labels)
        % Add rows and cells to the table and insert the images
        section = Section(['Phalange: ' params.phalange_labels(i_phalange)]);
        section.Numbered = false; % Remove section numbering
        
        tbl = Table();
        tbl.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
        for i = 1:2
            row = TableRow();
            for j = 1:2
                imgIndex = (j-1)*2 + i;
                img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections{imgIndex} '_butterfly_ph-' params.phalange_labels{i_phalange} '.jpg']));
                img.Style = {Width('8cm'), ScaleToFit};
                entry = TableEntry();
                append(entry, img);
                append(row, entry);
            end
            append(tbl, row);
        end
    
        add(section, tbl);
        add(chapter, section);
        add(chapter, PageBreak());
    end
    add(rpt, chapter);

    % Close the report
    close(rpt);
end
end