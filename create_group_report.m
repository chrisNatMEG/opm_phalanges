function create_group_report(save_path, subs, params, i_peak)
% Import necessary packages
import mlreportgen.report.*
import mlreportgen.dom.*

% Convert subNumber to a two-digit string
peak_label = ['_' params.peaks{i_peak}.label];

% Create the report name
reportName = ['group_report_' params.paradigm peak_label];

% Create a new report
rpt = Report(fullfile(save_path, reportName), 'pdf');

% Add a title page
titlepg = TitlePage;
titlepg.Title = ['Group results: ' params.paradigm];
add(rpt, titlepg);

% Add the table of contents
toc = TableOfContents();
add(rpt, toc);

% Define sections
sections = {'opm', 'squid', 'opmeeg', 'squideeg'};
sections2 = {'opm', 'squidmag', 'squidgrad'};

%% Params chapter
chapter = Chapter('Parameters');
processStruct(params, 'params', chapter);
add(rpt, chapter);

%% Sensor level 
chapter = Chapter('Sensor level');
chapter.Numbered = false; % Remove chapter numbering

dat = load(fullfile(save_path, ['group_sensor' peak_label]));
amp = dat.amp;
latency = dat.latency;
snr = dat.snr;
n_trl = dat.n_trl;

%% SNR_error section
section = Section('SNR_{error}');
section.Numbered = false; % Remove section numbering
dataCells = {};
dataCells{1} = snr.error_squidmag;
dataCells{2} = snr.error_opm;
rptTable = addGroupComparisonTables(dataCells, params.trigger_labels, {'squid', 'opm'}, '%.2f');
for i = 1:length(rptTable)
    add(section, rptTable{i});
    para = Paragraph(" ");
    add(section,para);
end
img = Image(fullfile(save_path,'figs',['SNR_error_box' peak_label '.jpg']));
img.Style = {Width('14cm'), ScaleToFit};
add(section, img);
add(chapter, section);

%% SNR_prestim section
section = Section('SNR_{prestim}');
section.Numbered = false; % Remove section numbering
dataCells = {};
dataCells{1} = snr.prestim_squidmag;
dataCells{2} = snr.prestim_opm;
rptTable = addGroupComparisonTables(dataCells, params.trigger_labels, {'squid', 'opm'}, '%.1f');
for i = 1:length(rptTable)
    add(section, rptTable{i});
    para = Paragraph(" ");
    add(section,para);
end
img = Image(fullfile(save_path,'figs',['SNR_prestim_box' peak_label '.jpg']));
img.Style = {Width('14cm'), ScaleToFit};
add(section, img);
add(chapter, section);

%% Amp section
section = Section('Amplitude (fT)');
section.Numbered = false; % Remove section numbering
dataCells = {};
dataCells{1} = amp.squidmag*1e15;
dataCells{2} = amp.opm*1e15;
rptTable = addGroupComparisonTables(dataCells, params.trigger_labels, {'squid', 'opm'}, '%.1f');
for i = 1:length(rptTable)
    add(section, rptTable{i});
    para = Paragraph(" ");
    add(section,para);
end
img = Image(fullfile(save_path,'figs',['Amplitude_meg_box' peak_label '.jpg']));
img.Style = {Width('14cm'), ScaleToFit};
add(section, img);
add(chapter, section);

%% Lat section
section = Section('Latency (ms)');
section.Numbered = false; % Remove section numbering
dataCells = {};
dataCells{1} = latency.squidmag*1e3;
dataCells{2} = latency.opm*1e3;
rptTable = addGroupComparisonTables(dataCells, params.trigger_labels, {'squid', 'opm'}, '%.1f');
for i = 1:length(rptTable)
    add(section, rptTable{i});
    para = Paragraph(" ");
    add(section,para);
end
img = Image(fullfile(save_path,'figs',['Latency_box' peak_label '.jpg']));
img.Style = {Width('14cm'), ScaleToFit};
add(section, img);
add(chapter, section);

%% Ntrl section
section = Section('N_{trls}');
section.Numbered = false; % Remove section numbering
dataCells = {};
dataCells{1} = n_trl.squidmag;
dataCells{2} = n_trl.opm;
rptTable = addGroupComparisonTables(dataCells, params.trigger_labels, {'squid', 'opm'}, '%.1f');
for i = 1:length(rptTable)
    add(section, rptTable{i});
    para = Paragraph(" ");
    add(section,para);
end
img = Image(fullfile(save_path,'figs',['Ntrl_MEG_box' peak_label '.jpg']));
img.Style = {Width('14cm'), ScaleToFit};
add(section, img);
add(chapter, section);

add(rpt, chapter);

%% Bad Channels chapter
% chapter = Chapter('Bad Channels');
% chapter.Numbered = false; % Remove chapter numbering
% for i_section = 1:length(sections)
%     filePath = fullfile(folderPath, ['sub_', subStr, '_' sections{i_section} '_badchs.mat']);
%     if isfile(filePath)
%         data = load(filePath);
%         section = Section(sections(i_section));
%         section.Numbered = false; % Remove section numbering
%         fields = fieldnames(data);
%         totalLength = 0;
%         for j = 1:length(fields)
%             varData = data.(fields{j});
%             if isempty(varData)
%                 varLength = 0;
%                 varList = '';
%             elseif iscell(varData)
%                 varLength = length(varData);
%                 varList = strjoin(varData, ', ');
%             else
%                 varLength = max(size(varData));
%                 varList = num2str(varData(:)');
%                 varList = regexprep(varList, '\s+', ', ');
%             end
%             totalLength = totalLength + varLength;
%             para = Paragraph([fields{j},' (n=', num2str(varLength), '): ', varList]);
%             add(section, para);
%         end
%         para = Paragraph(['n_{total}: ', num2str(totalLength)]);
%         add(section, para);
%         add(chapter, section);
%     end
% end
% add(rpt, chapter);

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
            if imgIndex <= length(sections2)
                img = Image(fullfile(save_path,'figs',[sections2{imgIndex} '_grndAvg_butterfly_trig-' params.trigger_labels{i_trigger} '.jpg']));
                img.Style = {Width('8cm'), ScaleToFit};
                entry = TableEntry();
                append(entry, img);
                append(row, entry);
            end
        end
        append(tbl, row);
    end
    add(section, tbl);
    add(chapter, section);
    add(chapter, PageBreak());
end
add(rpt, chapter);

%% Topographies
chapter = Chapter(['Topoplots']);
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
            if imgIndex <= length(sections2)
                imgIndex = (j-1)*2 + i;
                img = Image(fullfile(save_path,'figs',[sections2{imgIndex} peak_label '_grndAvg_topo_trig-' params.trigger_labels{i_trigger} '.jpg']));
                img.Style = {Width('8cm'), ScaleToFit};
                entry = TableEntry();
                append(entry, img);
                append(row, entry);
            end
        end
        append(tbl, row);
    end
    add(section, tbl);
    add(chapter, section);
    add(chapter, PageBreak());
end
add(rpt, chapter);

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