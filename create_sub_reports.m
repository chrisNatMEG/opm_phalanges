function create_sub_reports(subjectFolderPath, sub, params)
% Import necessary packages
import mlreportgen.report.*
import mlreportgen.dom.*

% Convert subNumber to a two-digit string
subStr = sprintf('%02d', sub);

% Create the report name
reportName = ['sub_', subStr, '_report'];

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
sections = {'opm', 'squid', 'opmeeg', 'squideeg'};
sections2 = {'opm', 'squidmag', 'squidgrad'};

%%
chapter = Chapter('Parameters');
processStruct(params, 'params', chapter);
add(rpt, chapter);

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
                varList = regexprep(varList, '\s+', ', ');
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
                varList = num2str(varData(:)');
                varList = regexprep(varList, '\s+', ', ');
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

        % Define the folder and the pattern
        pattern = ['sub_' subStr '_' sections{i_section} '_ica_rejected_comps*.jpg']; % Replace 'your_string' with the starting string
        files = dir(fullfile(subjectFolderPath, 'figs', pattern));

        if ~isempty(files)
            for i_file = 1:length(files)
                img = Image(fullfile(subjectFolderPath,'figs',files(i_file).name));
                img.Style = {ScaleToFit};
                add(section, img);
            end
        end
        add(chapter, section);
    end
end
add(rpt, chapter);


%% Butterflies
chapter = Chapter(['Butterfly plots']);
chapter.Numbered = false; % Remove chapter numbering
for i_trigger = 1:length(params.trigger_labels)
    % Add rows and cells to the table and insert the images
    section = Section(['Trigger: ' params.trigger_labels(i_trigger)]);
    section.Numbered = false; % Remove section numbering
    
    tbl = Table();
    tbl.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
    for i = 1:2
        row = TableRow();
        for j = 1:2
            imgIndex = (j-1)*2 + i;
            img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections{imgIndex} '_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']));
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

for i_peak = 1:length(params.peaks)
    %% chapter
    peak_label = params.peaks{i_peak}.label;
    chapter = Chapter(['Timelocked - ' peak_label]);
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
        data = load(fullfile(subjectFolderPath,['sub_' subStr '_' sections{i_section} '_' peak_label '.mat']));
        peak = data.peak;
        
        % Create the table with the required data
        num_triggers = length(params.trigger_labels);
        T = table('Size', [8, num_triggers], 'VariableTypes', repmat({'double'}, 1, num_triggers), 'VariableNames', params.trigger_labels);
        T.Properties.RowNames = {['peak_amplitude ' fieldUnit], ['max_amplitude ' fieldUnit], ['min_amplitude ' fieldUnit], 'peak_latency [ms]', 'SNR_prestim', 'SNR_stderr', ['std_prestim ' fieldUnit], ['stderr ' fieldUnit]};
        
        for i_trigger = 1:num_triggers
            data = peak{i_trigger};
            T{['peak_amplitude ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*data.peak_amplitude;
            T{['max_amplitude ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*data.max_amplitude;
            T{['min_amplitude ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*data.min_amplitude;
            T{'peak_latency [ms]', params.trigger_labels{i_trigger}} = 1e3*data.peak_latency;
            T{'SNR_prestim', params.trigger_labels{i_trigger}} = data.peak_amplitude / data.prestim_std;
            T{'SNR_stderr', params.trigger_labels{i_trigger}} = data.peak_amplitude / data.std_error;
            T{['std_prestim ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*data.prestim_std;
            T{['stderr ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*data.std_error;
        end
        
        % Convert the MATLAB table to a DOM table
        domTable = Table();
        domTable.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
        
        % Add header row
        headerRow = TableRow();
        append(headerRow, TableEntry(' '));
        for i_trigger = 1:num_triggers
            headerEntry = TableEntry(params.trigger_labels{i_trigger});
            headerEntry.Style = {HAlign('center'), Bold()};
            append(headerRow, headerEntry);
        end
        append(domTable, headerRow);
        
        formatString = {'%.1f','%.1f','%.1f','%.f','%.1f','%.2f','%.1f','%.1f'};

        % Add data rows
        for i_row = 1:height(T)
            row = TableRow();
            rowNameEntry = TableEntry(T.Properties.RowNames{i_row});
            rowNameEntry.Style = {Bold()};
            append(row, rowNameEntry);
            for j = 1:num_triggers
                entry = TableEntry(num2str(T{i_row, j},formatString{i_row}));
                entry.Style = {HAlign('right')}; 
                append(row, entry);
            end
            append(domTable, row);
        end
        add(section,domTable);
        add(chapter,section)
    end  
    
    %% Max channel plots
    for i_trigger = 1:length(params.trigger_labels)
        % Add rows and cells to the table and insert the images
        section = Section(['Peak channel - trigger: ' params.trigger_labels(i_trigger)]);
        section.Numbered = false; % Remove section numbering
        
        tbl = Table();
        tbl.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
        for i = 1:2
            row = TableRow();
            for j = 1:2
                imgIndex = (j-1)*2 + i;
                img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections{imgIndex} '_' peak_label '_evoked_peakchannel_trig-' params.trigger_labels{i_trigger} '.jpg']));
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

    %% Topo plots
    for i_trigger = 1:length(params.trigger_labels)
        % Add rows and cells to the table and insert the images
        section = Section(['Topoplot - phalange: ' params.trigger_labels(i_trigger)]);
        section.Numbered = false; % Remove section numbering
        
        tbl = Table();
        tbl.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
        for i = 1:2
            row = TableRow();
            for j = 1:2
                imgIndex = (j-1)*2 + i;
                img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections{imgIndex} '_' peak_label '_topo_trig-' params.trigger_labels{i_trigger} '.jpg']));
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

    %% Dipole plots
    if exist(fullfile(subjectFolderPath,[peak_label '_dipoles.mat']),'file')
        % Add rows and cells to the table and insert the images
        chapter = Chapter(['Dipoles - ' peak_label]);
        chapter.Numbered = false; % Remove section numbering
        
        for i_section = 1:length(sections2)
            section = Section(sections2(i_section));
            section.Numbered = false; % Remove section numbering

            img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections2{i_section} '_' peak_label '_dipfit_mri.jpg']));
            img.Style = {Width('14cm'), ScaleToFit};
            add(section, img);

            % Load the .mat file
            dipoles = load(fullfile(subjectFolderPath,[sections2{i_section} '_' peak_label '_dipoles.mat'])).dipoles;
            
            if params.numdipoles == 2
                % Create the table with the required data
                num_triggers = length(params.trigger_labels);
                T = table('Size', [4, num_triggers], 'VariableTypes', repmat({'double'}, 1, num_triggers), 'VariableNames', params.trigger_labels);
                T.Properties.RowNames = {['Dip_L mom [nAm]'], ['Dip_R mom [nAm]'], ['ResVar_L [%]'], ['ResVar_R [%]']};          
                for i_trigger = 1:num_triggers
                    data = dipoles{i_trigger};
                    T{['Dip_L mom [nAm]'], params.trigger_labels{i_trigger}} = 1e9*1e-4*max(vecnorm(data.dip.mom(1,:),2,1));
                    T{['Dip_R mom [nAm]'], params.trigger_labels{i_trigger}} = 1e9*1e-4*max(vecnorm(data.dip.mom(2,:),2,1));
                    T{['ResVar_L [%]'], params.trigger_labels{i_trigger}} = data.dip.rv(1)*100;
                    T{['ResVar_R [%]'], params.trigger_labels{i_trigger}} = data.dip.rv(1)*100;
                end
                
                % Convert the MATLAB table to a DOM table
                domTable = Table();
                domTable.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
                
                % Add header row
                headerRow = TableRow();
                append(headerRow, TableEntry(' '));
                for i_trigger = 1:num_triggers
                    headerEntry = TableEntry(params.trigger_labels{i_trigger});
                    headerEntry.Style = {HAlign('center'), Bold()};
                    append(headerRow, headerEntry);
                end
                append(domTable, headerRow);
                
                formatString = {'%.1f','%.1f','%.1f','%.1f'};
                % Add data rows
                for i_row = 1:height(T)
                    row = TableRow();
                    rowNameEntry = TableEntry(T.Properties.RowNames{i_row});
                    rowNameEntry.Style = {Bold()};
                    append(row, rowNameEntry);
                    for j = 1:num_triggers
                        entry = TableEntry(num2str(T{i_row, j},formatString{i_row}));
                        entry.Style = {HAlign('right')}; 
                        append(row, entry);
                    end
                    append(domTable, row);
                end

            else
                % Create the table with the required data
                num_triggers = length(params.trigger_labels);
                T = table('Size', [2, num_triggers], 'VariableTypes', repmat({'double'}, 1, num_triggers), 'VariableNames', params.trigger_labels);
                T.Properties.RowNames = {['Dip mom [nAm]'], ['ResVar [%]']};          
                for i_trigger = 1:num_triggers
                    data = dipoles{i_trigger};
                    T{['Dip mom [nAm]'], params.trigger_labels{i_trigger}} = fieldMultiplier*max(vecnorm(data.dip.mom(1,:),2,1));
                    T{['ResVar [%]'], params.trigger_labels{i_trigger}} = data.dip.rv(1)*100;
                end
                
                % Convert the MATLAB table to a DOM table
                domTable = Table();
                domTable.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
                
                % Add header row
                headerRow = TableRow();
                append(headerRow, TableEntry(' '));
                for i_trigger = 1:num_triggers
                    headerEntry = TableEntry(params.trigger_labels{i_trigger});
                    headerEntry.Style = {HAlign('center'), Bold()};
                    append(headerRow, headerEntry);
                end
                append(domTable, headerRow);
                
                formatString = {'%.1f','%.1f'};
    
                % Add data rows
                for i_row = 1:height(T)
                    row = TableRow();
                    rowNameEntry = TableEntry(T.Properties.RowNames{i_row});
                    rowNameEntry.Style = {Bold()};
                    append(row, rowNameEntry);
                    for j = 1:num_triggers
                        entry = TableEntry(num2str(T{i_row, j},formatString{i_row}));
                        entry.Style = {HAlign('right')}; 
                        append(row, entry);
                    end
                    append(domTable, row);
                end
            end
            add(section,domTable);
            add(chapter, section);
        end
        add(chapter, PageBreak());
        add(rpt,chapter);
    end

    %% MNE
    if exist(fullfile(subjectFolderPath,'opm_mne_peaks.mat'),'file')
        chapter = Chapter(['MNE - ' peak_label]);
        chapter.Numbered = false; % Remove section numbering
        
        for i_section = 1:length(sections2)
            section = Section(sections2(i_section));
            section.Numbered = false; % Remove section numbering

            tbl = Table();
            tbl.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
            for i = 1:ceil(length(params.trigger_labels)/2)
                row = TableRow();
                for j = 1:2
                    imgIndex = (j-1)*2 + i;
                    if imgIndex <= length(params.trigger_labels)
                        img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections2{i_section} '_' peak_label '_mne_trig-' params.trigger_labels{imgIndex} '.jpg']));
                        img.Style = {Width('8cm'), ScaleToFit};
                        entry = TableEntry();
                        append(entry, img);
                        append(row, entry);
                    end
                end
                append(tbl, row);
            end
            add(section, tbl);

            % Load the .mat file
            peaks = load(fullfile(subjectFolderPath,[sections2{i_section} '_mne_peaks.mat'])).peaks;
            
            if params.numdipoles == 2
                % Create the table with the required data
                num_triggers = length(params.trigger_labels);
                T = table('Size', [4, num_triggers], 'VariableTypes', repmat({'double'}, 1, num_triggers), 'VariableNames', params.trigger_labels);
                T.Properties.RowNames = {['FAHM_L [cm^2]'], ['FAHM_R [cm^2]'], ['Target region_L [%]'], ['Target region_R [%]']};          
                for i_trigger = 1:num_triggers
                    data = peaks{i_trigger};
                    T{['FAHM_L [cm^2]'], params.trigger_labels{i_trigger}} = data.fahm(1);
                    T{['FAHM_R [cm^2]'], params.trigger_labels{i_trigger}} = data.fahm(2);
                    T{['Target region_L [%]'], params.trigger_labels{i_trigger}} = data.target_region(1)*100;
                    T{['Target region_R [%]'], params.trigger_labels{i_trigger}} = data.target_region(2)*100;
                end
                
                % Convert the MATLAB table to a DOM table
                domTable = Table();
                domTable.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
                
                % Add header row
                headerRow = TableRow();
                append(headerRow, TableEntry(' '));
                for i_trigger = 1:num_triggers
                    headerEntry = TableEntry(params.trigger_labels{i_trigger});
                    headerEntry.Style = {HAlign('center'), Bold()};
                    append(headerRow, headerEntry);
                end
                append(domTable, headerRow);
                
                formatString = {'%.1f','%.1f','%.1f','%.1f'};
                % Add data rows
                for i_row = 1:height(T)
                    row = TableRow();
                    rowNameEntry = TableEntry(T.Properties.RowNames{i_row});
                    rowNameEntry.Style = {Bold()};
                    append(row, rowNameEntry);
                    for j = 1:num_triggers
                        entry = TableEntry(num2str(T{i_row, j},formatString{i_row}));
                        entry.Style = {HAlign('right')}; 
                        append(row, entry);
                    end
                    append(domTable, row);
                end

            else
                % Create the table with the required data
                num_triggers = length(params.trigger_labels);
                T = table('Size', [2, num_triggers], 'VariableTypes', repmat({'double'}, 1, num_triggers), 'VariableNames', params.trigger_labels);
                T.Properties.RowNames = {['FAHM [cm^2]'], ['Target region [%]']};          
                for i_trigger = 1:num_triggers
                    data = peaks{i_trigger};
                    T{['FAHM [cm^2]'], params.trigger_labels{i_trigger}} = data.fahm(1);
                    T{['Target region [%]'], params.trigger_labels{i_trigger}} = data.target_region(1)*100;
                end
                
                % Convert the MATLAB table to a DOM table
                domTable = Table();
                domTable.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
                
                % Add header row
                headerRow = TableRow();
                append(headerRow, TableEntry(' '));
                for i_trigger = 1:num_triggers
                    headerEntry = TableEntry(params.trigger_labels{i_trigger});
                    headerEntry.Style = {HAlign('center'), Bold()};
                    append(headerRow, headerEntry);
                end
                append(domTable, headerRow);
                
                formatString = {'%.1f','%.1f'};
    
                % Add data rows
                for i_row = 1:height(T)
                    row = TableRow();
                    rowNameEntry = TableEntry(T.Properties.RowNames{i_row});
                    rowNameEntry.Style = {Bold()};
                    append(row, rowNameEntry);
                    for j = 1:num_triggers
                        entry = TableEntry(num2str(T{i_row, j},formatString{i_row}));
                        entry.Style = {HAlign('right')}; 
                        append(row, entry);
                    end
                    append(domTable, row);
                end
            end
            add(section,domTable);
            add(chapter, section);
        end
        add(chapter, PageBreak());
        add(rpt, chapter);   
    end
end

%% Close the report
close(rpt);


%%% FUNCTIONS %%%%%%%%%%%%%%%%%%
function processStruct(s, parentName, chapter)
    localFields = fieldnames(s); % Make fields local
    for localI = 1:numel(localFields) % Make i local
        fieldName = localFields{localI};
        fieldValue = s.(fieldName);
        fullName = strcat(parentName, '.', fieldName);
        
        if isstruct(fieldValue)
            % If the field is a struct, process it recursively
            processStruct(fieldValue, fullName, chapter);
        elseif ismatrix(fieldValue) && isnumeric(fieldValue)
            % If the field is a 1-D array, convert it to a comma-separated string
            rowStrs = arrayfun(@(i) strjoin(arrayfun(@num2str, fieldValue(i,:), 'UniformOutput', false), ', '), 1:size(fieldValue,1), 'UniformOutput',false);
            valueStr = strjoin(rowStrs, '; ');
            append(chapter, Paragraph([fullName, ': ', valueStr]));
        elseif iscell(fieldValue) && all(cellfun(@ischar, fieldValue))
            % If the field is a cell array of strings, convert it to a comma-separated string
            valueStr = strjoin(fieldValue, ', ');
            append(chapter, Paragraph([fullName, ': ', valueStr]));
        elseif iscell(fieldValue) && all(cellfun(@isnumeric, fieldValue))
            % If the field is a cell array of strings, convert it to a comma-separated string
            valueStr = strjoin(cellfun(@num2str, fieldValue, 'UniformOutput', false), ', ');
            append(chapter, Paragraph([fullName, ': ', valueStr]));
        elseif isempty(fieldValue)
            % Otherwise, convert the value to a string and append it
            append(chapter, Paragraph([fullName, ': []']));
        elseif iscell(fieldValue) && ismatrix(fieldValue)
            for i_cell = 1:length(fieldValue)
                tmp = fieldValue{i_cell};
                processStruct(tmp,'tmp',chapter)
            end
        else
            % Otherwise, convert the value to a string and append it
            if size(fieldValue,1) > size(fieldValue,2)
                fieldValue = fieldValue';
            end
            valueStr = num2str(fieldValue);
            append(chapter, Paragraph([fullName, ': ', valueStr]));
        end
    end
end
end