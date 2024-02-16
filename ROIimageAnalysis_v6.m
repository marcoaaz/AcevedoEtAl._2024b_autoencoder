%% Root folder
clear clc
close all

scriptsFolder1 = 'E:\Alienware_March 22\scripts_Marco\updated MatLab scripts';
scriptsFolder2 = fullfile(scriptsFolder1, 'ROI');
scriptsFolder3 = fullfile(scriptsFolder1, 'WSI');
scriptsFolder4 = fullfile(scriptsFolder1, 'plots_miscellaneous');

addpath(scriptsFolder1);
addpath(scriptsFolder2);
addpath(scriptsFolder3);
addpath(scriptsFolder4);

%Gathering pre-registered images (phase map, PCA, DSA)

%Harzburgite
% dataFolder1 = 'E:\paper 3_datasets\harzburgite_synchrotron_christoph\tiff_HiRes\qupath_harzburgite\segmentation\11sep23_test4'; %semantic segmentation
% dataFolder2 = 'E:\paper 3_datasets\harzburgite_synchrotron_christoph\tiff_HiRes'; %X-ray stack
% dataFolder3 = 'E:\paper 3_datasets\harzburgite_synchrotron_christoph\tiff_HiRes\qupath_harzburgite'; %autoencoder image
% bad_list = {'sum', 'Ar', 'Cl', 'Flux0', 'Flux1', 'P', 'S', 'PbL', 'Si', 'Compton', 'elastic', 'Back'}; %Harzburgite 17-BSK-043

%Limestone
% dataFolder1 = 'E:\paper 3_datasets\conference paper_Goldschmidt\Export\registration\quPath_project\19Sep_RT_test3_original';
% dataFolder2 = 'E:\paper 3_datasets\conference paper_Goldschmidt\E3_66039_XFM scanning';
% expression1 = '\d*_(?<element>\w*)_\w*[.]tiff';
% bad_list = {'sum', 'Flux1', 'Al', 'Si', 'S', 'Cl', 'Br', 'Ar', 'Kr'}; %Limestone E3

%Gabbro
dataFolder1 = 'E:\paper 3_datasets\data_DMurphy\91714-81R5w-quartz-detail\tiff\quPath_gabbro\segmentation\11-sep_test3';
dataFolder2 = 'E:\paper 3_datasets\data_DMurphy\91714-81R5w-quartz-detail\tiff';
expression1 = '\d*-(?<element>\w*)[.]tiff';
bad_list = {'sum', 'Flux1', 'Flux0', 'Cl', 'Ar', 'Back'}; %Gabbro 81R-5-W






rootFolder = dataFolder1;
cd(rootFolder);

% DSAimage_rgb = imread(fullfile(dataFolder3, 'training_square_full_MedF_shortList4.tif'));
labelMap_rgb = imread(fullfile(dataFolder1, 'phasemap_target_RGB.tif'));
labelMap = imread(fullfile(dataFolder1, 'phasemap_target.tif')); %*_label.tif

[n_rows, n_cols] = size(labelMap);
img_size = [n_rows, n_cols];

%info
tablerank1 = readtable(fullfile(dataFolder1, 'species.xlsx'), 'Sheet', 'phaseMap');

minerals = tablerank1.Mineral; %PM_names = table2cell(tablerank(:, 2))'
triplet = [tablerank1.Triplet_1, tablerank1.Triplet_2, tablerank1.Triplet_3];
population = tablerank1.Pixels; 
n_masks = length(minerals);
density = repmat(2.7, [n_masks, 1]); 
%the density is required in Pie Chart by weight % (if known, import as a CSV file)

%Gathering quant map
suffix = '.tiff';
[fileNames, ~] = GetFileNames(dataFolder2, suffix); 


tokenNames = regexp(fileNames, expression1, 'tokens');
variableNames = string(tokenNames);
n_elements = length(variableNames);

quantMap = zeros(n_rows, n_cols, n_elements, 'double');
logMap = zeros(n_rows, n_cols, n_elements, 'double');
for i = 1:n_elements
    temp = imread(fullfile(dataFolder2, fileNames{i}));
    
    %medicine
    temp0 = temp;
    temp0(temp < 0) = 0; %clip zeros (for meaningful histograms and plots)
    
    temp1 = temp - min(temp, [], 'all'); %positive distribution (including negatives reduces noise)
    temp2 = log(temp1 + 1); %log-transform and prevent NaN (for meaningful dendrograms)
    
    %optional: percentile clipping
    
    %stack
    quantMap(:, :, i) = temp0; 
    logMap(:, :, i) = temp2;
end

%info
str_whos = whos('quantMap');
bytes_gb = (10^(-9))*str_whos.bytes; %compare w/ physical RAM memory
sprintf('quantMap consumes %.2f GB', bytes_gb)
memory 

variableNames %info print

%% Preparing data
close all

%Matrix --> Vector
dataLabels = variableNames;
indexVector = reshape(labelMap, [], 1);%maskPM vector for indexing

dataXYZ_quant = reshape(quantMap, [], n_elements); %unfolding data
dataXYZ_log = reshape(logMap, [], n_elements);

%Filtering out (from live feedback)

%Lines below N.B: for uncorrelated channels that have semi-quantitative data with
%artifacts (capping, etc.) or that represent mineral inclusions or grain
%boundary pixels, etc. (e.g.: S or Fe int-std in LA maps, , etc.). See
%below: It has to be defined for every mineral separately as a cell.
%default; 'S33', 'Cu65', 'Zn66', 'Pb208', 'As75'

%Option 1: using an undesired list
% bad_list = {'elastic', 'Compton', 'Al', 'sum', 'Flux1'}; 

% bad_list = {'sum', 'Flux1', 'Flux0', 'Cl', 'Ar', 'Back'}; %Gabbro 81R-5-W
% bad_list = {'sum', 'Ar', 'Cl', 'Flux0', 'Flux1', 'P', 'S', 'PbL', 'Si', 'Compton', 'elastic', 'Back'}; %Harzburgite 17-BSK-043
bad_list = {'sum', 'Flux1', 'Al', 'Si', 'S', 'Cl', 'Br', 'Ar', 'Kr'}; %Limestone E3

bad_idx = contains(dataLabels, bad_list);
good_idx = ~bad_idx; %option 1: using bad_list

%Option 2: using good_list
%good_list = {'K', 'Ca'};
% good_list = {'Ni', 'Cr', 'Ca', 'Fe'};
% good_idx = contains(dataLabels, good_list);

dataLabels2 = dataLabels(good_idx);
dataXYZ_quant = dataXYZ_quant(:, good_idx);
dataXYZ_log = dataXYZ_log(:, good_idx);
dataXYZ2_log = normalize(dataXYZ_log, 1, 'range');

%Print mineral list (for next section)
choicesTable_full = cell2table(minerals, 'VariableNames', {'Minerals'}, ...
    'RowNames', string(1:n_masks));
disp(choicesTable_full);

%message
text1 = sprintf('The elements will be: %s', strjoin(string(dataLabels2), ', '));
text1

%% Optional: Pre-check images
%Noisy layers should be included in bad_list

str_ele = 'Back' %element (choose manually)
str_idx = strcmp(dataLabels2, str_ele); %if none, change element

%To see outliers
temp_quantMap = dataXYZ_quant(:, str_idx); 
temp_img = reshape(temp_quantMap, img_size(1), img_size(2));

%percentile (too see structure), optional
percentOut = 0.005;
P = prctile(temp_img, [percentOut, 100 - percentOut], "all");
temp_img = rescale(temp_img, 0, 1, 'InputMin', P(1), 'InputMax', P(2));           

temp_img2 = uint8(255*rescale(temp_img, 0, 1)); %issue with NaNs
colormap_temp = colormap("jet");

figure 
imshow(temp_img2, colormap_temp)
title(str_ele)

figure, 
histogram(temp_quantMap)
title(str_ele)



%% ROI Tool parameters (desired list)

ref_img = labelMap_rgb;
% ref_img = DSAimage_rgb;
% ref_img = [labelMap_rgb; DSAimage_rgb];

%remember: pairs_registered = {QM_ref, PM_RGB_reg, BSE_reg, TM_ref_reg};

%defining ROIs
loup_x = n_cols/2; %default location
loup_y = n_rows/2;

sideLength = 6; 
zoomLength = 100;
box_pos = [...
    ceil(loup_x) - sideLength/2, ceil(loup_y) - sideLength/2, ...
    sideLength, sideLength];
zoom_pos = [...
    ceil(loup_x) - zoomLength/2, ceil(loup_y) - zoomLength/2, ...
    zoomLength, zoomLength];
%[xmin, ymin, width, height]; upper left corner (<6^2 pixels is fast).

radius = 100; %probe

%plot settings
th_bluePlaneTH = 60; %for sparsity plot
th_alpha = 0.5;
plotOptions = [1]; %pie chart (area=1 & weight=2)
hist_xAxisLabel = '';
histogramScale = "linear";
dim_biplot = 3; %for contribution plot & biplot
dim_dendrogram = 4;
th_dendrogram = 0.5; %pct of maximum dendrogram height

%Working arrays (keeping all is computationally expensive)
choices = [1, 2, 3, 4];
maskLabel = minerals(choices); %minerals(choices, 1);

choicesTable = cell2table(maskLabel, 'VariableNames', {'Minerals'}, ...
    'RowNames', string(choices));
disp(choicesTable);
choice = input('Select mineral number (for contribution plot):\n'); %for ROIdimGraphs

%binscatter plot axis (of the two, none should be empty)
x1 = 'Ca';  %Fe
x2 = 'Fe';  %Cr
pctOut = 0.05; %fot knowing binscatter plot axis limits

%messages
text1 = sprintf('The elements will be: %s', strjoin(string(dataLabels2), ', '));
text2 = sprintf('The minerals will be: %s', strjoin(string(maskLabel), ', '));
text3 = sprintf('The 1D histograms will plot: %s', string(maskLabel(choices == choice)));
text4 = sprintf('The biscatter plots will display: X= %s vs Y= %s', x1, x2);

{text1; text2; text3; text4}

%ROIClicked, ROIMoved, MovingROI

%% Probe 1: 
% linear data: null matrix, pie chart, binscatter, histograms
%log-transformed data: scree plot, contribution plot, biplot, dendrogram

close all

figure

imshow(ref_img); %use reference image
ax = gca;

hCircle = images.roi.Ellipse('Center', [loup_x loup_y], ...
    'Semiaxes', [radius radius], ...
    'Parent', ax, ...
    'Color', 'r');

%pre-allocate
[hCircle] = dataROI(hCircle, img_size, dataXYZ_quant, dataXYZ2_log, ...
    indexVector, minerals, choices, dataLabels2, "");

addlistener(hCircle, 'MovingROI',...
    @(varargin)dataROI(hCircle, img_size, dataXYZ_quant, dataXYZ2_log, ...
    indexVector, minerals, choices, dataLabels2, ""));

% %informative
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIbarGraph3D(hCircle, th_bluePlaneTH, th_alpha)); %change faceAlpha > 0 for the blue plane
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIpieChart(hCircle, minerals, density, triplet, plotOptions)); 

%pre-selected list: 2D histograms
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIbinscatterGraph(hCircle, dataLabels2, minerals, choices, ...
    x1, x2, pctOut)); %pre-selected minerals

%choose one: 1D histogram
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIhistogramsGraph(hCircle, choices, dataLabels2, 3*pctOut, ...
    hist_xAxisLabel, histogramScale));
% 
% %Probe 2: selected minerals 
%the input to dataROI is log-transformed
addlistener(hCircle, 'MovingROI',...
    @(varargin)dataROI(hCircle, img_size, dataXYZ_quant, dataXYZ2_log, ...
    indexVector, minerals, choices, dataLabels2, ""));

%choose one: scree/contribution plot
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIdimGraph(hCircle, choice, dataLabels2, 1:dim_biplot));

%pre-selected list: dendrogram
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIbiplotGraph(hCircle, choices, dataLabels2, 1:dim_biplot)); %requires all components of pca: 'Economy'
addlistener(hCircle, 'ROIMoved',...
    @(varargin)ROIdendrogramPCA(hCircle, choices, dataLabels2, 1:dim_dendrogram, 'ward', th_dendrogram)); 

%% Probe 3 (optional): pixel quality control (QA/QC)
%linear data: quality control probe

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

addlistener(hProbe, 'MovingROI',...
    @(varargin)dataROI(hProbe, img_size, dataXYZ_quant, dataXYZ2_log, ...
    indexVector, minerals, choices, dataLabels2, ""));
%2nd half of function (nullMatrix): dataLabels, minerals and lastColumns
%{'BSE_bf1'} if present ; ("" if empty) necessary

%Quality control probe %xyzModality
addlistener(hProbe, 'ROIMoved',...
    @(varargin)ROIpixelQC(hProbe, dataXYZ_quant, dataLabels2)); 

%% Optional: 

close all

hFig = figure;
hFig.WindowState = 'maximized';

t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% First plot
ax1 = nexttile;
imshow(labelMap_rgb)

% Second plot
ax2 = nexttile;
imshow(DSAimage_rgb)

linkaxes([ax1 ax2], 'xy')

hTitle = sgtitle('Synchronized view');
hTitle.FontSize = 12;