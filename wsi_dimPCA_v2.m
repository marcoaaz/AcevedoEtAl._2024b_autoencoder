%Dimensionality reduction: first script
%Author: Marco Andres, ACEVEDO ZAMORA
%Citation: https://doi.org/10.1016/j.chemgeo.2024.121997

%Creation: 13-11-2023, Marco Acevedo
%Update: 29-05-2024, Marco Acevedo

%%
clear 
clc

%Required libraries
scriptPath1 = 'E:\Alienware_March 22\current work\00-new code May_22\';
scriptPath2 = 'E:\Alienware_March 22\scripts_Marco\updated MatLab scripts';
addpath(fullfile(scriptPath1))
addpath(fullfile(scriptPath2))
addpath(fullfile(scriptPath1, '\rayTracing'))
addpath(fullfile(scriptPath1, '\dimReduction\'))
addpath(fullfile(scriptPath2, '\external_package\'))

%*********** User input ***************

workingDir = 'D:\Chris_Collaboration\2024_Ioan_Purdys_reward\XFM\153874_P2\tiff\hybrid\recoloured_trial1_pctOut0.1';
cd(workingDir)

%Input List
%note 1: adjust based on retrospective insight (the first should include everyting)

%153874_P2 --> general
undesiredList = {'AgL', 'Ar', 'CdL', 'Cl', ...
    'count_rate', 'dead_time', 'dwell', 'elastic', 'flux', ...
    'Flux0', 'Flux1', 'K', 'Kr', 'P', 'pileup', 'raw_flux', ...
    'Rb', 'S', 'Si', 'sum'};

desiredList = {'Chlorite_Fe', 'Chlorite_Compton', 'Chlorite_Back', ...
    'Chlorite_Ca', 'Chlorite_count_rate', 'Chlorite_K', 'Chlorite_Ti', ...
    'Chlorite_Y', 'Chlorite_Zr', 'Au_AuL'};
% desiredList = {'Fe', 'Compton', 'Transmission', 'Back', ...
%     'Ca', 'count_rate', 'K', 'Ti', ...
%     'Y', 'Zr', 'AuL'};

%Imaging parameters: 
percentOut_pcaImage = 1; %default = 3;

%Define saving names
imageVer = '19-Jun-24'; %files destDir
denomination1 = 'trial1'; %name of trial run


%*************** Script *******************

stats_table = readtable(fullfile(workingDir, "descriptiveStats.csv"));
experiment_names = stats_table.experiment;
element_names = stats_table.element;
layer_names = strcat(experiment_names, '_', element_names);
%%
folder1 = 'original_pyramid';
folder2 = 'linear_pyramid';
original_path = fullfile(workingDir, strcat(folder1, '_files'), "/0/");
recoloured_path = fullfile(workingDir, strcat(folder2, '_files'), "/0/");
recoloured_xml_path = fullfile(workingDir, strcat(folder2, '.dzi'));

%Gathering montage info

[content] = xml2struct(recoloured_xml_path); %*.dzi format
img_format = content.Image.Attributes.Format;
imageHeight = str2double(content.Image.Size.Attributes.Height);
imageWidth = str2double(content.Image.Size.Attributes.Width);
tileWidth = str2double(content.Image.Attributes.TileSize);
tileHeight = tileWidth;

tileNames = GetFileNames(recoloured_path, img_format); 
n_tiles = length(tileNames); %within montage
info_st = imfinfo(fullfile(recoloured_path, tileNames{1}));
temp = info_st.BitsPerSample;
n_layers = length(temp); %SamplesPerPixel
bitDepth = temp(1);

%sorting
name_pattern = '\d*_\d*.tif'; %default (edit manually)
fileNames2 = regexp(tileNames, name_pattern, 'match');
fileNames3 = string(fileNames2(~any(cellfun('isempty', fileNames2), 1)));
q0 = regexp(fileNames3, '\d*', 'match'); %sort, respect 'down & right' order
q1 = str2double(cat(1, q0{:}));
[q1_index, ii] = sortrows(q1, [1, 2]); %edit manually (ascending +)

tileNames_sorted = tileNames(ii); %input file names
%number of tiles in each montage row/col
n_rows_montage = max(unique(q1_index(:, 2))) + 1; 
n_cols_montage = max(unique(q1_index(:, 1))) + 1;

%Generate folders
tag_combo = strcat(imageVer, '_', denomination1);
destDir = fullfile(workingDir, strcat('pca_results/', tag_combo));
destDir1 = fullfile(destDir, 'pca_tiles');
mkdir(destDir1);

%Input filter
n_names = length(layer_names); %parsed with regex in 'tilingAndStacking_v3.py'

%finding undesired
idx_out = ones([1, n_names]); %1= not to include
if exist("desiredList", 'var')    
    
    for i = 1:n_names %finding desired  
        
        for j = 1:length(desiredList)
            cond1 = strfind(layer_names{i}, desiredList{j});            
            if cond1 == 1                                     
                idx_out(i) = 0;
                break
                %warning: 'S' might be confused with 'Si' or 'Sr'
                %idx_out 
            else                    
                idx_out(i) = 1;                
            end
        end
    end        
else    
    
    for i = 1:n_names %finding undesired    
        for j = 1:length(undesiredList)
            cond2 = isempty(strfind(layer_names{i}, undesiredList{j}));
            if cond2
                idx_out(i) = 0;
                %warning: 'S' might be confused with 'Si' or 'Sr'
            else
                idx_out(i) = 1;
                break
            end
        end
    end    
end
idx_in = ~idx_out;
desiredList = layer_names(idx_in);


%%
%PCA the incremental way
disp('Running incremental PCA')

n_layers_input = sum(idx_in);
for i = 1:n_tiles    
    fprintf('tile %d/%d (x=%d_y=%d)  \n', ...
        i, n_tiles, q1_index(i, 1), q1_index(i, 2))

    fileName_input = fullfile(recoloured_path, tileNames_sorted{i});
    stack_input = imread(fileName_input); %int8
    stack_subset = stack_input(:, :, idx_in);
    n_layers_input2 = size(stack_subset, 3);

    mtx_long = double( transpose(reshape(stack_subset, ...
        [], n_layers_input2)) );    

    if i == 1 %conventional
        n_samples_seen = size(mtx_long, 2);    
        mu = mean(mtx_long, 2);
        [U, S, ~] = svd(mtx_long - mu, 'econ');
        n_components = []; %to consider in calculation    

        continue
    end
    %loop
    [U, S, mu, n_samples_seen] = incrementalPCA...
           (mtx_long, U, S, mu, n_samples_seen, n_components);    
end 
fprintf('%s completed.\n', 'pca')

%Storing data
vect_max = max(S, [], 1);
vect_sum = vect_max/sum(S, 'all')*100;
energy_thresh = sum(vect_sum(1:3)); %explained

pca_montage.U_full = U; %PC's (columns)
pca_montage.S = S; %SVD diagonal (columns, sorted)
pca_montage.mu = mu; %mean, single column
pca_montage.n_samples_seen = n_samples_seen;
pca_montage.energy_thresh = energy_thresh;

%info
idx_search = find(cumsum(vect_sum) > 95, 1);
sprintf('3PC explain: %.1f pct', energy_thresh)
sprintf('95pct of variability explained by %.0f PCs \n out of %.0f selected variables', ...
    idx_search, n_layers_input)

%Saving process metadata (reproducibility)
s.workingDir = workingDir;
s.destDir = destDir1;
s.desiredList = desiredList; %from idx_in
s.undesiredList = undesiredList;
s.idx_in = idx_in; %referenced to 'descriptiveStats.csv' rows
s.imageHeight = imageHeight;
s.imageWidth = imageWidth;
s.tileHeight = tileHeight;
s.tileWidth = tileWidth;
s.percentOut_pcaImage = percentOut_pcaImage;
s.pca_montage = pca_montage;

save(fullfile(destDir, 'montage_pcaInfo.mat'), 's', '-mat', '-v7.3', "-nocompression");

%Apply PCA and generate reconstruction image
disp('Applying PCA eigen-vectors')
U_sub = U(:, 1:3); %for RGB

M = 2; %number of workers, total=8 with 16 threads (logical processors)
%Note: monitor that there is not RAM overload to reduce M.
parfor (i = 1:n_tiles, M)
    img_stats = [] %clearing (prevents runtime error)
    
    fileName_input = fullfile(recoloured_path, tileNames_sorted{i});
    stack_input = imread(fileName_input); %int8
    stack_subset = stack_input(:, :, idx_in);
    
    %properties
    temp_height = size(stack_subset, 1);
    temp_width = size(stack_subset, 2);
    n_layers_input2 = size(stack_subset, 3);
    
    X = double( transpose(reshape(stack_subset, ...
        [], n_layers_input2)) ); %mtx_long    
    X_demean = X - mu;
    score_3pc = X_demean'*U_sub;
    
    %Reconstructing
    pc1 = reshape(score_3pc(:, 1), temp_height, temp_width); %R
    pc2 = reshape(score_3pc(:, 2), temp_height, temp_width); %G
    pc3 = reshape(score_3pc(:, 3), temp_height, temp_width); %B
    
    img_stats = single(cat(3, pc1, pc2, pc3));
            
    %Configure file saving
    fullFileName = fullfile(destDir1, strcat(tileNames_sorted{i}));

    t = Tiff(fullFileName, 'w');
    tagstruct = [];
    tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 3;
    tagstruct.SampleFormat = 3;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';        
    tagstruct.ImageLength = temp_height;
    tagstruct.ImageWidth = temp_width; 
        
    setTag(t, tagstruct)
    write(t, img_stats);
    close(t);  
end

%Building mosaic (for visualisation)
img_montage = zeros(imageHeight, imageWidth, 3, 'single');
for i = 1:n_tiles

    temp_tile_name = fullfile(destDir1, tileNames_sorted{i});
    temp_tile = imread(temp_tile_name);
    
    from_row = 1 + q1_index(i, 2)*tileHeight;
    to_row = min(from_row + tileHeight - 1, imageHeight);
    from_col = 1 + q1_index(i, 1)*tileWidth;
    to_col = min(from_col + tileWidth - 1, imageWidth);

    % disp([from_row, to_row, from_col, to_col])
    img_montage(from_row:to_row, from_col:to_col, :) = temp_tile;
    
end

%re-scaling
pc1 = img_montage(:, :, 1);
pc2 = img_montage(:, :, 2);
pc3 = img_montage(:, :, 3);
P1 = prctile(pc1, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
P2 = prctile(pc2, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
P3 = prctile(pc3, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
pc1_rs = rescale(pc1, 0, 255, 'InputMin', P1(1), 'InputMax', P1(2));
pc2_rs = rescale(pc2, 0, 255, 'InputMin', P2(1), 'InputMax', P2(2));
pc3_rs = rescale(pc3, 0, 255, 'InputMin', P3(1), 'InputMax', P3(2));

img_double = cat(3, pc1_rs, pc2_rs, pc3_rs);
pca_rgb = uint8(img_double);

%fixing image
%pca_rgb2 = flipud(pca_rgb); %due to pyvips convention
% pca_rgb = imnlmfilt(pca_rgb); %denoising

%saving image
file_name1 = fullfile(destDir, strcat('incrementalPCA_rgb_pctOutput', num2str(percentOut_pcaImage),'.tif'));
imwrite(pca_rgb, file_name1);

close all
figure, 
imshow(pca_rgb)
