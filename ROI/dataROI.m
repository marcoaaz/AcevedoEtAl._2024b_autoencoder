
function [roiHandle] = dataROI(roiHandle, img_size, dataXYZ, dataXYZ_log, ...
    maskPM, minerals, choices, dataLabels, lastColumns)

%Input:

%roiHandle = drawn ROI
%img_size = size of image reference (equal to image stack)
%dataXYZ/_log = unfolded stack
%maskPM = phase map as vector;
% minerals = string with names of all original phases ranked following
% label numbers in qupathPhaseMap_v7.m

%% XYZ generator 

n_masks = length(minerals);

logicalMask = createMask(roiHandle, img_size(1), img_size(2));
roiDataXYZ = dataXYZ(logicalMask(:), :); %for quantification
roiDataXYZ_log = dataXYZ_log(logicalMask(:), :); %for multi-variate statistics
roiMineralMask = maskPM(logicalMask(:)); 

pixelPopulations = zeros(1, n_masks);
for i = 1:n_masks
    pixelPopulations(i) = sum(roiMineralMask == i);
end

roiHandle.UserData.logicalMask = logicalMask(:);
roiHandle.UserData.roiDataXYZ = roiDataXYZ;
roiHandle.UserData.roiMineralMask = roiMineralMask;
roiHandle.UserData.pixelPopulations = pixelPopulations;

%% Null matrix
%matrix: minerals x variables

temp4 = cell(1, n_masks);
for i = 1:n_masks
   temp1 = roiDataXYZ(roiMineralMask == i, :);
   temp2 = 100*sum(temp1 == 0, 1)/sum(roiMineralMask == i); %Zero pixels 
   temp3 = 100*sum(isnan(temp1), 1)/sum(roiMineralMask == i); %Nan pixels
   
   temp4{i} = array2table(temp2 + temp3); %for step below   
end
nullPixels = vertcat(temp4{:}); %concatenating table rows vertically
nullPixels.Properties.VariableNames = dataLabels;

if ~(lastColumns == "") %rearranging columns
    newOrder = {setdiff(dataLabels, lastColumns), lastColumns}; %placing last
    nullPixels = nullPixels(:, [newOrder{:}]);    
end
nullPixels = addvars(nullPixels, minerals, 'Before', 1); %placing first

roiHandle.UserData.nullPixels = nullPixels;

%% Statistical analysis 

% %Relevant Pearson coefficients
% [lkup] = rankedPearson(maskXYZ, dataLabels, 0.6, 0.05); %expensive computation
% roiHandle.UserData.lkup = lkup;

%HCA of PCA (selected mineral masks)
maskLabel = minerals(choices, 1);
pixelNumber = cell(1, length(choices));
coefs = cell(1, length(choices));
score = cell(1, length(choices));
S = cell(1, length(choices));
k = 0;
for i = choices
    k = k + 1;
    
    logical = roiMineralMask == i;
    pixelNumber{k} = sum(logical);
       
    maskXYZ = roiDataXYZ_log(logical, :);
    
    %PCA
    [coefs{k}, score{k}, S{k}] = pca(maskXYZ, 'Economy', true, 'Rows', 'all');
    % [coefs{k}, score{k}, S{k}] = pca(maskXYZ, 'Economy', false, 'Rows',
    % 'pairwise'); %uses Eigenvalue decomposition not SVD loadings
    % , for the n (observations)-by-p(variables)
    % coefficient is p-by-p for each PC in descending order of variance
end

roiHandle.UserData.minerals = minerals;
roiHandle.UserData.choices = choices;
roiHandle.UserData.maskLabel = maskLabel;
roiHandle.UserData.pixelNumber = pixelNumber;
roiHandle.UserData.coefs = coefs;
roiHandle.UserData.score = score;
roiHandle.UserData.S = S;

end
