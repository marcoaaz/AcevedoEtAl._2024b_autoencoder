%ROIimageAnalysis_v7_wsi.m

%It allows loading a quantitative map image stack for analysis with a
%region of interest tool (ROI Tool) that is hoovered on top of a reference
%image (phase map, PCA, DSA map) for interrogating different minerals
%(separately) and their compositional variations and multi-variate analysis.

%Updated: 15-Jun-24, Marco Acevedo
%Citation: https://www.sciencedirect.com/science/article/pii/S0009254124000779?via%3Dihub


%% Root folder
clear 
clc
close all

scriptsFolder1 = 'E:\Alienware_March 22\scripts_Marco\updated MatLab scripts';
scriptsFolder2 = fullfile(scriptsFolder1, 'ROI');
scriptsFolder3 = fullfile(scriptsFolder1, 'ROI_wsi');
scriptsFolder4 = fullfile(scriptsFolder1, 'plots_miscellaneous/');
scriptsFolder5 = fullfile(scriptsFolder1, 'external_package/');
addpath(scriptsFolder1);
addpath(scriptsFolder2);
addpath(scriptsFolder3);
addpath(scriptsFolder4)
addpath(scriptsFolder5)

%User input 
%Note: ROI Tool cannot do multi-variate statistics with 1 element only

%Gabbro 81R-5-W
dataFolder1 = 'E:\paper 3_datasets\data_DMurphy\91714-81R5w-quartz-detail\tiff\quPath_gabbro\segmentation\11-sep_test3';
dataFolder2 = 'E:\paper 3_datasets\data_DMurphy\91714-81R5w-quartz-detail\tiff\recoloured_trial1_pctOut0.1';
expression1 = '\d*-(?<element>\w*)[.]tiff';
suffix = '.tiff';
bad_list = {'sum', 'Flux1', 'Flux0', 'Cl', 'Ar', 'Back'}; %gathered from retrospective insight
good_list = {'Cr', 'Ca', 'Fe', 'P', 'Si', 'Mn', 'K'}; %only used when defined

%probes settings
sideLength = 8; 
zoomLength = 100;
radius = 100; %probe

%plot settings
th_bluePlaneTH = 60; %for sparsity plot
th_alpha = 0.5;
plotOptions = [1]; %pie chart (area=1 & weight=2)
hist_xAxisLabel = '';
histogramScale = "linear";
dim_biplot = 3; %for contribution plot & biplot
dim_dendrogram = 4; %>4 elements selected
th_dendrogram = 0.5; %pct of maximum dendrogram height
%if the biplot elements are too unrelated the dendrograms blank (height = 1).

%binscatter plot axis (they are in the 'good_list')
x1 = 'Ca';  %Fe
x2 = 'Si';  %Cr
pctOut = 0.05; %for knowing binscatter plot axis limits

%Minerals to focus interrogation (speeds up interaction)
choices = [2, 3]; %default [1, 2, 3, 4]

%Script

%Follows 'tilingAndStacking_v3.py' convention
original_xml_path = fullfile(dataFolder2, 'recoloured_pyramid.dzi');
stats_path = fullfile(dataFolder2, 'descriptiveStats.csv');
original_path = fullfile(dataFolder2, 'original_pyramid_files', '0');
recoloured_path = fullfile(dataFolder2, 'recoloured_pyramid_files', '0');

%Pre-aligned segmentation data ('qupathPhaseMap_v8.m' convention)
labelMap = imread(fullfile(dataFolder1, 'phasemap_target.tif')); %*_label.tif
labelMap_rgb = imread(fullfile(dataFolder1, 'phasemap_target_RGB.tif'));

ref_img = labelMap_rgb; %DSAimage_rgb , [labelMap_rgb; DSAimage_rgb];
img_size = size(ref_img);

%metadata
tablerank1 = readtable(fullfile(dataFolder1, 'species.xlsx'), 'Sheet', 'phaseMap');
minerals = tablerank1.Mineral;
n_masks = length(minerals);
triplet = [tablerank1.Triplet_1, tablerank1.Triplet_2, tablerank1.Triplet_3];
population = tablerank1.Pixels; 
density = repmat(2.7, [n_masks, 1]); %for pie chart, import as a CSV file

%message 1
choicesTable_full = cell2table(minerals, 'VariableNames', {'Minerals'}, ...
    'RowNames', string(1:n_masks));
disp(choicesTable_full);

%Quantitative maps info 
stats_table = readtable(stats_path);
element_names = stats_table.element;
variableNames = string(element_names);
n_elements = length(variableNames);

%Filtering layers out (from live feedback)
bad_idx = contains(variableNames, bad_list);
if ~exist('good_list', 'var')
    good_idx = ~bad_idx; %option 1: using bad_list
else
    good_idx = contains(variableNames, good_list);
end
variableNames2 = variableNames(good_idx);

%Montage info
[content] = xml2struct(original_xml_path); %*.dzi format
img_format = content.Image.Attributes.Format;
imageWidth = str2double(content.Image.Size.Attributes.Width);
imageHeight = str2double(content.Image.Size.Attributes.Height);
tileWidth = str2double(content.Image.Attributes.TileSize);
tileHeight = tileWidth;

%The ROI Tool displays a smaller front-end image for greater speed
scale_factor = imageHeight/img_size(1); 

tileNames = GetFileNames(original_path, img_format); 
n_tiles = length(tileNames); %within montage
info_st = imfinfo(fullfile(original_path, tileNames{1}));
temp = info_st.BitsPerSample;
n_layers = length(temp); %SamplesPerPixel
bitDepth = temp(1);

%tile name sorting
name_pattern = '\d*_\d*.tif'; %default (edit manually)
fileNames2 = regexp(tileNames, name_pattern, 'match');
fileNames3 = string(fileNames2(~any(cellfun('isempty', fileNames2), 1)));
q0 = regexp(fileNames3, '\d*', 'match'); %sort, respect 'down & right' order
q1 = str2double(cat(1, q0{:}));
[q1_index, ii] = sortrows(q1, [1, 2]); %edit manually (ascending +)

%number of tiles in each montage row/col
n_rows_montage = max(unique(q1_index(:, 2))) + 1; 
n_cols_montage = max(unique(q1_index(:, 1))) + 1;
tileNames_sorted = tileNames(ii); %input file names

%Tile corners (in pixels)
corners_array = zeros([n_tiles, 4], "double");
k = 0;
for j = 1:n_cols_montage %down-right order
    for i = 1:n_rows_montage
        k = k + 1;
        tl_row = 1 + (i-1)*tileHeight;
        tl_col = 1 + (j-1)*tileWidth;
        br_row = min(tl_row + tileHeight - 1, imageHeight);
        br_col = min(tl_col + tileWidth - 1, imageWidth);

        corners_array(k, :) = [tl_row, tl_col, br_row, br_col];
    end
end
n_rows_vec = corners_array(:, 3) - corners_array(:, 1) + 1;
n_cols_vec = corners_array(:, 4) - corners_array(:, 2) + 1; 
corners_array_WH = [corners_array, n_rows_vec, n_cols_vec];
corners_array_WH_ds = floor(corners_array_WH/scale_factor); %factors of 2 work

%Defining ROI Tool
loup_x = img_size(2)/2; %default location
loup_y = img_size(1)/2;
box_pos = [...
    ceil(loup_x) - sideLength/2, ceil(loup_y) - sideLength/2, ...
    sideLength, sideLength];
zoom_pos = [...
    ceil(loup_x) - zoomLength/2, ceil(loup_y) - zoomLength/2, ...
    zoomLength, zoomLength];
%[xmin, ymin, width, height]; upper left corner (<6^2 pixels is fast).

%Choices for analysis
maskLabel = minerals(choices); %minerals(choices, 1);
choicesTable = cell2table(maskLabel, 'VariableNames', {'Minerals'}, ...
    'RowNames', string(choices));
disp(choicesTable);
choice = input('Select mineral number (for contribution plot):\n'); %for ROIdimGraphs

%message 3
text0 = 'During the analysis';
text1 = sprintf('The elements will be: %s', strjoin(string(variableNames2), ', '));
text2 = sprintf('The minerals will be: %s', strjoin(string(maskLabel), ', '));
text3 = sprintf('The 1-D histograms will use: %s', string(maskLabel(choices == choice)));
text4 = sprintf('The biscatter plots will use: X= %s vs Y= %s', x1, x2);
{text0; text1; text2; text3; text4}

%structure input for callback functions
roiTool_metadata.original_path = original_path; %image pyramids
roiTool_metadata.recoloured_path = recoloured_path;
roiTool_metadata.minerals = minerals; %segmentation metadata
roiTool_metadata.n_masks = n_masks;
roiTool_metadata.triplet = triplet;
roiTool_metadata.population = population;
roiTool_metadata.density = density;
roiTool_metadata.variableNames2 = variableNames2; %quantification metadata
roiTool_metadata.choices = choices; %index of minerals to focus
roiTool_metadata.good_idx = good_idx; %index of elements to focus
roiTool_metadata.img_size = img_size; %size of reference image
roiTool_metadata.imageWidth = imageWidth; %size of montage
roiTool_metadata.imageHeight = imageHeight;
roiTool_metadata.tileWidth = tileWidth; %size of montage tiles
roiTool_metadata.tileHeight = tileHeight;
roiTool_metadata.n_tiles = n_tiles;
roiTool_metadata.n_layers = n_layers;
roiTool_metadata.n_rows_montage = n_rows_montage;
roiTool_metadata.n_cols_montage = n_cols_montage;
roiTool_metadata.tileNames_sorted = tileNames_sorted;
roiTool_metadata.corners_array_WH = corners_array_WH;
roiTool_metadata.corners_array_WH_ds = corners_array_WH_ds;
roiTool_metadata.scale_factor = scale_factor;

%% ROI Tool
%actions allowed: ROIClicked, ROIMoved, MovingROI

tile_str = strcat('Tile', {' '}, string(1:n_tiles));

%Probe 1
% linear data: null matrix, pie chart, binscatter, histograms
% log-transformed data: scree plot, contribution plot, biplot, dendrogram

clc
close all

figure
ax = gca;

imshow(ref_img); %use reference image

%draw tiles for reference (optional)
for s = 1:n_tiles
    pos_temp = [corners_array(s, 2) - 0.5, corners_array(s, 1) - 0.5,...
        tileWidth, tileHeight];
    
    images.roi.Rectangle('Position', pos_temp, 'Parent', ax, ...
        'Color', 'white', 'EdgeAlpha', 0.3, 'FaceAlpha', 0);    
end

text( (corners_array(:, 2)+ corners_array(:, 4))/2, ...
    (corners_array(:, 1)+ corners_array(:, 3))/2, ...
    tile_str, 'Color', 'black', 'FontSize', 30, 'FontWeight','bold', ...
    'HorizontalAlignment','center');

%Draw
hCircle = images.roi.Ellipse('Center', [loup_x loup_y], ...
    'Semiaxes', [radius radius], ...
    'Parent', ax, ...
    'Color', 'yellow');

%pre-allocate
[hCircle] = dataROI_wsi(hCircle, roiTool_metadata, labelMap, "");

addlistener(hCircle, 'MovingROI',...
    @(varargin)dataROI_wsi(hCircle, roiTool_metadata, labelMap, ""));

%informative
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIbarGraph3D(hCircle, th_bluePlaneTH, th_alpha)); %change faceAlpha > 0 for the blue plane
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIpieChart2(hCircle, roiTool_metadata, plotOptions)); 

%pre-selected list: 2D histograms
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIbinscatterGraph2(hCircle, roiTool_metadata, ...
    x1, x2, pctOut)); %pre-selected minerals

%choose one: 1D histogram
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIhistogramsGraph2(hCircle, roiTool_metadata, 3*pctOut, ...
    hist_xAxisLabel, histogramScale));

%  Probe 2: selected minerals

%2- or 3-D biplot
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIbiplotGraph2(hCircle, roiTool_metadata, 1:dim_biplot)); 

%4-D dendrogram (pre-selected list)
%chosen one: scree/contribution plot
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIdimGraph(hCircle, choice, variableNames2, 1:dim_dendrogram));

addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIdendrogramPCA2(hCircle, roiTool_metadata, ...
    1:dim_dendrogram, 'ward', th_dendrogram)); 

%% Probe 3 (optional): pixel quality control (QA/QC)
%applied on linear data: quality control probe

clc

figure

imshow(ref_img, 'InitialMagnification', 550); %use reference image
xlim([zoom_pos(1), zoom_pos(1) + zoom_pos(3)])
ylim([zoom_pos(2), zoom_pos(2) + zoom_pos(4)])

ax = gca; %alternative: fig = gcf; ax = fig.CurrentAxes;

hProbe = drawrectangle(...
    ax,...
    'Position', box_pos,...
    'StripeColor', 'r',...
    'InteractionsAllowed', 'translate'); %images.roi.Rectangle

%pre-allocate
[hProbe] = dataROI_wsi(hProbe, roiTool_metadata, labelMap, "");

addlistener(hProbe, 'MovingROI',...
    @(varargin)dataROI_wsi(hProbe, roiTool_metadata, labelMap, ""));

%Quality control probe %xyzModality
addlistener(hProbe, 'ROIMoved',...
    @(varargin)ROIpixelQC2(hProbe, roiTool_metadata)); 

