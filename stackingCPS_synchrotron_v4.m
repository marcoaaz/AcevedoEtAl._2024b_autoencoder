%Dimensionality reduction SCRIPT 1
%Author: Marco Andres, ACEVEDO ZAMORA
%Last update: 13-11-2023
%Citation: 

clear 
clc

%Libraries
scriptPath1 = 'E:\Alienware_March 22\current work\00-new code May_22\';
scriptPath2 = 'E:\Alienware_March 22\scripts_Marco\updated MatLab scripts';
addpath(fullfile(scriptPath1, '\dimReduction\'))
addpath(fullfile(scriptPath2))

%****************** User input ***************
%Importing info
inputType = 'edx'; %depending on the stack, change between 'edx'/'xfm'
workingDir = 'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\Balz_Rock report';
fNames = GetFileNames(workingDir, '.png'); 
%Note: % xfm='.tiff'; % edx= '.png'

undesiredList = {'BSE'}; %adjust based on retrospective insight

% %Limestone E3
% undesiredList = {'Ar', 'Br', 'Ge', 'Kr', 'PbL', 'Se', 'As', 'Zr', 'HfL', ...
% 'elastic', 'sum', 'Si', 'Flux1'}; 
% undesiredList = strcat(undesiredList, '_NoNaNs.tiff'); %note: needs to match exactly 

%Imaging parameters: 
percentOut = 0.1; %default= 0.5 for XFM; 0.1 for log-transformed XFM
filterSize = 1; %default= [3, 3] neighbourhood for EDX; 5 for XFM
portion = 100; %area pct for PCA (default = 10%)
percentOut_pcaImage = 1; %default = 3;

%Saving names
imageVer = '6-Dec'; %files destDir
denomination1 = 'test1'; %name of trial run

tic;

%*************** Script routine *******************
cd(workingDir)

denomination2 = sprintf('%dx%d', filterSize, filterSize);
denomination3 = strcat('outInput', num2str(percentOut));
denomination4 = strcat('outPCA', num2str(percentOut_pcaImage));

% gathering info
info_st = imfinfo(fullfile(workingDir, fNames{1}));
imageHeight = info_st.Height;
imageWidth = info_st.Width;
stackType = {'edx', 'xfm'};
idxType = find(strcmp(stackType, inputType));
if idxType == 1 %optional: EDX accessory phases might dissappear
    percentOut = 0;
end

%Calculate image size to use
imageDim_max = 10000; %default = 10000, if larger than original preserves original 
%Note: working with a 10K cap prevents looping with >4GB stacks and
%overloading the RAM during the Python workflow

[val_max, idx_max] = max([imageHeight, imageWidth]);
if (val_max > imageDim_max) && (idx_max == 1)    
    imageDim_use = [imageDim_max, NaN];

elseif (val_max > imageDim_max) && (idx_max == 2)
    imageDim_use = [NaN, imageDim_max];

elseif (val_max <= imageDim_max)
    imageDim_use = [imageHeight, imageWidth];

else 
    imageDim_use = [imageDim_max, imageDim_max];
end

%Destination directories
destDir = fullfile(workingDir, strcat(inputType, '_', imageVer)); %'synchrotron_stack_manualFilter'
rescaled_dir = fullfile(destDir, denomination1);
mkdir(destDir);
mkdir(rescaled_dir);

%Script convention: Intermediate file naming
name_str = strcat(imageVer, ...
    '_', denomination1, ...
    '_', denomination2, ...
    '_', denomination3);

fileNameRescaled = strcat(name_str, '.tif'); %'_rs_denoised_short.tif';
fileName_pcaDescrip = strcat(name_str, '_', denomination4, '_pcaInfo', '.mat');% 'pca_struct.mat'
fileNamePCA = strcat(name_str, '_', denomination4, '_pcaRGB', '.tif'); %'denoised_pca_short.tif';
fileNameForAutoencoder = strcat(name_str, '_imageStack', '.mat'); %'stack_temp_ext_MedF_short.mat';

%Filtering out undersired layers
n_names = length(fNames);
n_out = length(undesiredList);
idx_out = true([1, n_names]);
for i = 1:n_names
    for j = 1:n_out

        if isempty(strfind(fNames{i}, undesiredList{j}))
            idx_out(i) = 0;
            %Note of caution: 'S' might be confused with 'Si' or 'Sr'
        else
            idx_out(i) = 1;
            break
        end
    end
end
fNames1 = fNames(~idx_out);
n_types = length(fNames1);
fNames1'

%Image stacking
stack_temp = [];
for i = 1:n_types
    imageName = fNames1{i};    

    switch idxType
        case 1
            %for SEM-EDX (16-bit TIMA exports)
            temp1 = imread(fullfile(workingDir, imageName)); %edx
            temp1 = temp1(:, :, 1);
        case 2
            %for Synchrotron
            temp1 = imread(fullfile(workingDir, imageName)); %synchrotron  

            %Optional: Medicine 
            % temp1(temp1 < 0) = 0; %clip zeros: optional (improves PC pct explaining variables)
            
            %Test
            temp1 = temp1 - min(temp1, [], 'all'); %positive distribution (reduces noise)
            temp1 = log(temp1 + 1); %log-transform and prevent NaN (dendrograms)
    end    
    temp1 = imresize(temp1, imageDim_use); %Optional: Resizing image 
  
    %Rescaling (Option 1): For Synchrotron
    P = prctile(temp1, [percentOut, 100 - percentOut], "all");
    temp2 = rescale(temp1, 0, 255, 'InputMin', P(1), 'InputMax', P(2));           

    %Optional denoising (not for low resolution 1000x500; only high 3501x1501)
    temp3 = medfilt2(temp2, [filterSize, filterSize]); %admits double   
    stack_temp = cat(3, stack_temp, temp3);%Stacking     

    %Saving (for obs.)
    newName = strrep(imageName, '.tiff', fileNameRescaled);
    imwrite(uint8(temp3), jet(256), fullfile(rescaled_dir, newName), 'compression', 'none')
end

%Save (for Autoencoder script)
destFile = fullfile(destDir, fileNameForAutoencoder);
save(destFile, 'stack_temp', '-mat', '-v7.3') %mat for binary 

t.import= toc;

%% Principal component analysis in RGB

[n_rows, n_cols, n_channels] = size(stack_temp);
XY_size = n_rows*n_cols;
image_XYZ3 = reshape(stack_temp, XY_size, n_channels, []); %[] position changes how data fold
X_full = double(image_XYZ3); 
%Note 1: use 'double' rather than 'single' if the following messate
%appears: Warning: Columns of X are linearly dependent to within machine
%precision. Using only the first 1 components to compute TSQUARED. 

% X_full = normalize(X_full, 1, 'range'); %optional
%Note 2: image layer normalisation is optional since it might worsen noise
%in the PCA image. 

% Note 3: There is no problem in not normalising since the later DSA script
% in Python and 'ROIimageAnalysis_v6.m' dendrograms automatically do
% normalisation

%PCA calculation
s_seed = rng;
rng(s_seed)
%random sample: 
% 12 sec for 4000x6500 px; 
% 30sec 17Kx9Kx21

idx = randperm(XY_size); 
from = round((100 - portion)*XY_size/100) + 1; %alternative: use a window of 5Kx5K
X = X_full(idx(from:end), :); 
[coeff, score, latent, ~, explained, mu] = pca(X); %centers variables automatically

t.pca= toc;

%PCA Projection (~3min for 10Kx15K)
score_3pc = (X_full - mu)*coeff(:, 1:3); %obtaining PCs scores

%Reconstructing: t = score1*coeff1' + repmat(mu1,13,1)
pc1 = reshape(score_3pc(:, 1), n_rows, n_cols); %R
pc2 = reshape(score_3pc(:, 2), n_rows, n_cols); %G
pc3 = reshape(score_3pc(:, 3), n_rows, n_cols); %B

%re-scaling
P1 = prctile(pc1, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
P2 = prctile(pc2, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
P3 = prctile(pc3, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
pc1_rs = rescale(pc1, 0, 255, 'InputMin', P1(1), 'InputMax', P1(2));
pc2_rs = rescale(pc2, 0, 255, 'InputMin', P2(1), 'InputMax', P2(2));
pc3_rs = rescale(pc3, 0, 255, 'InputMin', P3(1), 'InputMax', P3(2));

img_double = cat(3, pc1_rs, pc2_rs, pc3_rs);
pca_rgb = uint8(img_double);
clear pc1 pc2 pc3

t.reprojecting= toc;

%optional
% pca_rgb = imnlmfilt(pca_rgb); 
t.denoising = toc;

%Quick info:
close all
figure, 
imshow(pca_rgb)

%Command Window info
sprintf('3PC explain: %.1f pct', sum(explained(1:3)))
idx_search = find(cumsum(explained) > 95, 1);
sprintf('95pct of variability explained by %.0f PCs', idx_search)

%% Saving process information

%Info for reproducing the image analysis outputs
s.workingDir = workingDir;
s.allFiles = fNames;
s.usedFiles = fNames1;
s.destDir = destDir;
s.rescaled_dir = rescaled_dir;
s.imageDim_max = imageDim_max;
s.imageHeight = imageHeight;
s.imageWidth = imageWidth;
s.percentOut = percentOut;
s.percentOut_pcaImage = percentOut_pcaImage;
s.filterSize = filterSize;
s.inputType = inputType;
s.portion = portion;
s.seed = s_seed;
s.coeff = coeff;
s.explained = explained;
s.mu = mu;
save(fullfile(destDir, fileName_pcaDescrip), '-struct', 's');

%saving image
imwrite(pca_rgb, fullfile(destDir, strcat(imageVer, '_', fileNamePCA)));


% %% Optional 1: Check stack layers as recoloured images
% 
% test = temp1(:, :, 1);
% 
% figure
% imshow(uint16(temp1)) %after importing
% figure
% imshow(uint8(temp2)) %after rescaling
% figure
% imshow(uint8(temp3)) %after median filter
% 
% sel = 6;%2, 5, end
% 
% figure
% imshow(stack_temp(:, :, sel), 'Colormap', jet(256))
% colorbar 
% title(fNames1{sel})
% 
% figure
% histogram(temp2)
% title(fNames1{sel})
% 
% %% Optional 2: Pre-check before PCA image production
% close all
% 
% sel = 1;
% fNames1{sel}
% 
% figure
% imshow(uint8(stack_temp(:, :, sel)), jet(256))
% 
% figure,
% histogram(X_full(:, sel))


