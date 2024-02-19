%% Root directories
clear 
clc

workingDir = 'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\GeoPIXE meeting';
cd(workingDir)

phasemapFile = 'phasemap_target_RGB.tif';
labelFile = 'phasemap_label.tif';
infoFile = 'species.xlsx';
%only one to change
% phasemapDir = 'E:\paper 3_datasets\data_DMurphy\91714-81R5w-quartz-detail\tiff\quPath_gabbro\segmentation\11sep23_trial1'; %gabbro
phasemapDir = 'E:\paper 3_datasets\harzburgite_synchrotron_christoph\tiff_HiRes\qupath_harzburgite\segmentation\11sep23_test4'; %harzburgite
% phasemapDir = 'E:\paper 3_datasets\conference paper_Goldschmidt\Export\registration\quPath_project\19Sep_RT_test3'; %limestone

datasetDescription = 'harzburgite';
currentFile = 'test2';

imgRGB = imread(fullfile(phasemapDir, phasemapFile));
imgLabel = imread(fullfile(phasemapDir, labelFile));
infoTable = readtable(fullfile(phasemapDir, infoFile), 'Sheet', 'phaseMap');

originalLabels = infoTable.NewLabel;
minerals = infoTable.Mineral;
minerals %printed info
dim_ref = size(imgRGB);
n_rows = dim_ref(1);
n_cols = dim_ref(2); %check w/ debugg section if resized

%% Reading GeoPIXE region

fileName = 'Marco-region0.csv'; %original 'boilerplate' region to modify
text_2 = fileread(fileName);

expression2 = '\n[a-zA-Z_]+[,]';
c = regexp(text_2, expression2, 'start');
d = regexp(text_2, expression2, 'end');

n_characters = length(text_2);
n_items = length(c);

header_str = string(n_items);
data_cell = cell(1, n_items);
idx_mtx = zeros(n_items, 4, 'int64');
for i = 1:n_items
    %headers 
    header_start = c(i) + 1;
    header_end = d(i) - 1;    
    data_start = d(i) + 1;
    try
        data_end = c(i + 1) - 1;
    catch
        data_end = n_characters - 1;
    end

    temp_header = text_2(header_start:header_end);           
    temp_str = strsplit( convertCharsToStrings(text_2(data_start:data_end)) , ',');
    
    header_str(i) = temp_header;
    data_cell{i} = temp_str;
    idx_mtx(i, :) = [header_start, header_end, data_start, data_end];
end
% %debugg
% dim = str2double(data_cell{4});
% n_rows = dim(2); %y
% n_cols = dim(1); %x
% data_cell
% header_str'

%% Drawing Target pixels

close all

hFig = figure;
imshow(imgRGB)
ax = gca;
% hROI = drawcircle("Center", uint16([dim_ref(2), dim_ref(1)]/2), 'Radius', 150, ...
%     'parent', ax);

%full area
x_tl = 0;
y_tl = 0;
box_w = n_cols;
box_h = n_rows;

hROI = drawrectangle("Position", [x_tl, y_tl, box_w, box_h], 'Parent', ax);

[hROI] = maskWithinROI(hROI, [dim_ref(1), dim_ref(2)], imgLabel(:), minerals); %initial update

addlistener(hROI, 'MovingROI',...
    @(varargin)maskWithinROI(hROI, [dim_ref(1), dim_ref(2)], imgLabel(:), minerals));

%% Saving GeoPIXE files

destDir = fullfile(workingDir, datasetDescription);
mkdir(destDir)

pixelPopulations = hROI.UserData.pixelPopulations;
indList = hROI.UserData.indList;

idx_contained = (pixelPopulations ~= 0);
indList2 = indList(idx_contained);
minerals_contained = minerals(idx_contained);

n_contained = sum(idx_contained);
for m = 1:n_contained

    sel_str = minerals_contained{m};
    filename = fullfile(destDir, strcat(currentFile, '_', sel_str, '.csv'));  % better to use fullfile(path,name) 
    
    %Image dim
    char_dim = char(strjoin(string([n_cols, n_rows]), ','));

    %Q-vector
    idx_target = indList2{m}; %col-major order
    [col, row] = ind2sub([n_rows, n_cols], idx_target); %conversion
    idx_target2 = sub2ind([n_cols, n_rows], row, col);
    idx_target3 = sort(idx_target2, 'ascend');
    char_ROI = char(strjoin(string(idx_target3), ','));
    
    %introduced changes
    change_list = {'Image', 'Note', 'Q'};
    change_items = {char_dim, char(sel_str), char_ROI};
    
    writeGeopixeRegion(text_2, header_str, data_cell, idx_mtx, change_list, change_items, filename)

end

%optional
% logical_mask = false(n_rows, n_cols);
% logical_mask(idx_target) = 1;
% figure, imshow(logical_mask)
% imwrite(logical_mask, 'target_logicalMask.tif');

%% Optional: Visualise GeoPIXE region (2D image, not RGB)
%figure flipped upside down in GeoPIXE

x = str2double(data_cell{23}); %centroids
y = str2double(data_cell{24}); 
ind_img = str2double(data_cell{25});%row-major idx
[col, row] = ind2sub([n_cols, n_rows], ind_img); %conversion

img_rescaled = imresize(imgRGB, [n_rows, NaN]);
img_rescaled(row, col) = 255;

figure
imshow(img_rescaled)
hold on
plot(x, y, '*')
hold off
%%

function writeGeopixeRegion(text_2, header_str, data_cell, idx_mtx, change_list, change_items, filename)
% description: Writing GeoPIXE region

%Input

%a : GeoPIXE text file (comment section)
%header_str : items available (medata can be re-updated within GeoPIXE)
%data_cell : values of those items
%idx_mtx : char idx ranges of headers and data (table: n_headers x 4)
%change_list : headers that we want to change
%change_items : values of those headers
%filename : looped file name (for each mineral)

%Output

%saved file to import in https://asci.synchrotron.org.au/

%%

expression1 = '#';
a = regexp(text_2, expression1, "end"); %for rewriting file

n_changes = length(change_list);
n_items = length(header_str);

%update ranges according to new data
idx_mtx_accum = idx_mtx;
data_cell_updated = data_cell;

for ii = 1:n_changes
    change_idx = find(ismember(header_str, change_list{ii}));
    temp_char_idx = idx_mtx_accum(change_idx, :); %row to change
    replacement_char = change_items{ii};
    
    original_length = 1 + temp_char_idx(4) - temp_char_idx(3);
    dif_length = original_length - length(replacement_char);    
    
    idx_mtx_accum(change_idx:end, :) = idx_mtx_accum(change_idx:end, :) - dif_length;
    idx_mtx_accum(change_idx, 3) = temp_char_idx(:, 3); %data start preserved

    data_cell_updated{change_idx} = change_items{ii};
end

%comments
line1 = a(end);
char_1 = [text_2(1:line1), newline];

char_2 = char();
for jj = 1:n_items
    temp_data_str = data_cell_updated{jj};
    temp_data_char = char(strjoin(string(temp_data_str), ','));

    char_2 = [char_2, char(header_str{jj}), ',', temp_data_char, newline];
end

char_full = [char_1, char_2];

fid = fopen(filename, 'w');    % open file for writing (overwrite if necessary)
fprintf(fid,'%s', char_full);          % Write the char array, interpret newline as new line
fclose(fid);                  % Close the file (important)
% open(filename) 

end

function [roiHandle] = maskWithinROI(roiHandle, img_size, maskPM, minerals)

%Input:

%roiHandle = drawn ROI
%img_size = size of image reference (equal to image stack)
%maskPM = phase map as vector;
% minerals = string with names of all original phases ranked following
% label numbers in qupathPhaseMap_v7.m

n_masks = length(minerals);
logicalMask = createMask(roiHandle, img_size(1), img_size(2));
roiMineralMask = maskPM(logicalMask(:)); 

pixelPopulations = zeros(1, n_masks);
for i = 1:n_masks
    pixelPopulations(i) = sum(roiMineralMask == i);
end

%indices
indList = cell(1, n_masks);
for j = 1:n_masks
    temp_logical = find( (maskPM == j) & logicalMask(:) );
    indList{j} = temp_logical;
end

roiHandle.UserData.logicalMask = logicalMask(:);
roiHandle.UserData.roiMineralMask = roiMineralMask;
roiHandle.UserData.pixelPopulations = pixelPopulations;
roiHandle.UserData.indList = indList;

% assignin('base', "indList", indList)
end
