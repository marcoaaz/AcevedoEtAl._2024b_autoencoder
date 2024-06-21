%Dimensionality reduction: first script
%Author: Marco Andres, ACEVEDO ZAMORA
%Citation: https://doi.org/10.1016/j.chemgeo.2024.121997

%Creation: 13-11-2023, Marco Acevedo
%Update: 15-06-2024, Marco Acevedo

clear 
clc

%****************** User input ***************

%Required libraries
scriptPath1 = 'E:\Alienware_March 22\current work\00-new code May_22\';
scriptPath2 = 'E:\Alienware_March 22\scripts_Marco\updated MatLab scripts';

%Importing info
dirList = {
    'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\Appendix_C_Tornado_Maps\Lady_Annie_Tornado\MS_LAN028\elements';
    'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\Appendix_C_Tornado_Maps\Lady_Annie_Tornado\MS_LAN032\elements';
    'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\Appendix_C_Tornado_Maps\Pegmont_Tornado\MS_PEG025B\elements';
    'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\Appendix_C_Tornado_Maps\Pegmont_Tornado\MS_PEG028\elements';
    'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\Appendix_C_Tornado_Maps\Pegmont_Tornado\MS_PEG039\element';
    'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\Appendix_C_Tornado_Maps\Pegmont_Tornado\MS_PEG050\elements';
    };
list_n = length(dirList);

undesiredList = {''};

%undesiredList = strcat(undesiredList, '_NoNaNs.tiff');

%note 1: adjust based on retrospective insight (the first should include everyting)
%note 2: 'undesiredList' needs to match exactly (check for original name typos)

%Imaging parameters: 
inputType = 'xfm'; %depending on the stack, change between 'edx'/'xfm'
img_format = '.png'; %xfm='.tiff'; % edx= '.png'
percentOut = 0.1; %default= 0.5 for XFM; 0.1 for log-transformed XFM
filterSize = 3; %default= [3, 3] neighbourhood for EDX; 5 for XFM
portion = 100; %area pct for PCA (default = 10%)
percentOut_pcaImage = 1; %default = 3;

%Define saving names
imageVer = '15-Jun-24'; %files destDir
denomination1 = 'trial1'; %name of trial run

%*************** Script *******************

tic;
addpath(fullfile(scriptPath1))
addpath(fullfile(scriptPath2))
addpath(fullfile(scriptPath1, '\rayTracing'))
addpath(fullfile(scriptPath1, '\dimReduction\'))

%Loop implementation
for k = 1:list_n
    workingDir = dirList{k};

    %Destination directories
    destDir = fullfile(workingDir, strcat(inputType, '_', imageVer)); %'synchrotron_stack_manualFilter'
    rescaled_dir = fullfile(destDir, denomination1);
    mkdir(destDir);
    mkdir(rescaled_dir);
    cd(workingDir)
    
    %Script convention: Intermediate file naming
    denomination2 = sprintf('%dx%d', filterSize, filterSize);
    denomination3 = strcat('outInput', num2str(percentOut));
    denomination4 = strcat('outPCA', num2str(percentOut_pcaImage));
    
    name_str = strcat(imageVer, ...
        '_', denomination1, ...
        '_', denomination2, ...
        '_', denomination3);
    
    fileNameRescaled = strcat(name_str, '.tif'); %'_rs_denoised_short.tif';
    fileName_pcaDescrip = strcat(name_str, '_', denomination4, '_pcaInfo', '.mat');% 'pca_struct.mat'
    fileNamePCA = strcat(name_str, '_', denomination4, '_pcaRGB', '.tif'); %'denoised_pca_short.tif';
    fileNameForAutoencoder = strcat(name_str, '_imageStack', '.mat'); %'stack_temp_ext_MedF_short.mat';
    
    % gathering info
    fNames = GetFileNames(workingDir, img_format); 
    info_st = imfinfo(fullfile(workingDir, fNames{1}));
    imageHeight = info_st.Height;
    imageWidth = info_st.Width;
    
    %Filtering out undersired layers
    n_names = length(fNames);
    n_out = length(undesiredList);
    idx_out = true([1, n_names]);
    for i = 1:n_names
        disp(i)
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
    
    %settings
    stackType = {'edx', 'xfm'};
    idxType = find(strcmp(stackType, inputType));
    if idxType == 1 %optional, since EDX accessory phases might dissappear
        percentOut = 0;
    end
    
    %Capping image size (to avoid RAM overload)
    imageDim_max = 8000; %default = 8000, %if larger than original preserves original 
    %Note 1: ~8K cap (in test) prevents PCA step crash
    %Note 2: 10K cap prevents loop output >4GB stacks that crash the autoencoder (next script)
    
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
    
    %Image stacking
    stack_temp = [];
    for i = 1:n_types
        imageName = fNames1{i};    
        temp = double(imread(fullfile(workingDir, imageName)));
        
        if size(temp, 3) > 1
            temp0 = temp(:, :, 1);
        end
    
        switch idxType
            case 1 %for SEM-EDX (16-bit TIMA exports)            
                temp1 = temp0;
    
            case 2 %for Synchrotron              
                %medicine (optional)
                % temp0(temp0 < 0) = 0; 
                %clip zeros (optional): improves principal component 
                % pct explaining variables
                
                %log-transform and prevent NaN (for dendrograms)
                temp1_a = temp0 - min(temp0, [], 'all'); %positive distribution (reduces noise)
                temp1 = log(temp1_a + 1); 
    
        end    
    
        temp1_b = imresize(temp1, imageDim_use); %resizing (optional)
      
        %Rescaling (option for Synchrotron)
        P = prctile(temp1_b, [percentOut, 100 - percentOut], "all");
        temp2 = rescale(temp1_b, 0, 255, 'InputMin', P(1), 'InputMax', P(2));           
    
        %Optional denoising (not for low resolution 1000x500; only high 3501x1501)
        temp3 = medfilt2(temp2, [filterSize, filterSize]); %admits double   
        
        %Saving full image (for retrospective feedback)
        newName = strrep(imageName, '.tiff', fileNameRescaled);
        imwrite(uint8(temp3), jet(256), ...
            fullfile(rescaled_dir, newName), 'compression', 'none')
    
        stack_temp = cat(3, stack_temp, temp3); %stacking     
    
        
    
        sprintf('%.f completed = %s', i, imageName);
    end
    clear temp0 temp1 temp1_a temp1_b temp2 temp3 
    
    t.importing= toc;
    
    tic;
    %Save (for Autoencoder script)
    destFile = fullfile(destDir, fileNameForAutoencoder);
    save(destFile, 'stack_temp', '-mat', '-v7.3') %mat for binary 
    t.saving = toc;
    
    tic;
    %Principal component analysis in RGB        
    
    [n_rows, n_cols, n_channels] = size(stack_temp);
    XY_size = n_rows*n_cols;
    image_XYZ3 = reshape(stack_temp, XY_size, n_channels, []); %[] position changes how data folds
    X_full = double(image_XYZ3); 
    %Note 1: use 'double' rather than 'single' if the following message appears: 
    %Warning: Columns of X are linearly dependent to within machine
    %precision. Using only the first 1 components to compute TSQUARED. 
    
    %Normalization (optional)
    % X_full = normalize(X_full, 1, 'range');
    %Note 2: image layer normalisation is optional since it might worsen noise
    %in the PCA image. 
    % Note 3: There is no problem in not normalising since the later DSA script
    % in Python and 'ROIimageAnalysis_v6.m' dendrograms do automatic normalisation
    
    s_seed = 1;
    rng(s_seed) %random number generator
    
    idx = randperm(XY_size); 
    from = round((100 - portion)*XY_size/100) + 1; %alternative: use a window of 5Kx5K
    
    X = X_full(idx(from:end), :); 
    [coeff, score, latent, ~, explained, mu] = pca(X); %centers variables automatically
    
    t.pca= toc;
    %PCA calculation: 12 sec for 4000x6500 px; 30sec for 17Kx9Kx21
    
    tic;
    %PCA Projection (~3min for 10Kx15K)
    score_3pc = (X_full - mu)*coeff(:, 1:3); %obtaining PCs scores
    
    %Reconstructing: t = score1*coeff1' + repmat(mu1,13,1)
    pc1 = reshape(score_3pc(:, 1), n_rows, n_cols); %R
    pc2 = reshape(score_3pc(:, 2), n_rows, n_cols); %G
    pc3 = reshape(score_3pc(:, 3), n_rows, n_cols); %B
    
    %Re-scaling
    P1 = prctile(pc1, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
    P2 = prctile(pc2, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
    P3 = prctile(pc3, [percentOut_pcaImage, 100 - percentOut_pcaImage], "all");
    pc1_rs = rescale(pc1, 0, 255, 'InputMin', P1(1), 'InputMax', P1(2));
    pc2_rs = rescale(pc2, 0, 255, 'InputMin', P2(1), 'InputMax', P2(2));
    pc3_rs = rescale(pc3, 0, 255, 'InputMin', P3(1), 'InputMax', P3(2));
    
    img_double = cat(3, pc1_rs, pc2_rs, pc3_rs);
    pca_rgb = uint8(img_double);
    clear pc1 pc2 pc3 pc1_rs pc2_rs pc3_rs
    
    t.reprojecting= toc;
    
    tic;
    %optional
    % pca_rgb = imnlmfilt(pca_rgb); 
    t.denoising = toc;
    t    

    %Command Window info    
    sprintf('During PCA, 3 PC are explaining: %.1f pct', sum(explained(1:3)))
    idx_search = find(cumsum(explained) > 95, 1);
    sprintf('95pct of variability explained by %.0f PCs', idx_search)
    
    %Saving process metadata (reproducibility)
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
    
end

%% Quality check
close all
figure, 
imshow(pca_rgb)
