
function [roiHandle] = dataROI_wsi(roiHandle, roiTool_metadata, ...
    labelMap, lastColumns)

%Input
%roiHandle = drawn ROI

img_size = roiTool_metadata.img_size; %'ref_img' (important: it fits in RAM)
n_tiles = roiTool_metadata.n_tiles;
minerals = roiTool_metadata.minerals;
choices = roiTool_metadata.choices;
dataLabels = roiTool_metadata.variableNames2;

n_masks = length(minerals);
logicalMask = createMask(roiHandle, img_size(1), img_size(2)); %'ref_img'

%FoV
try 
    roi_condition_prev = roiHandle.UserData.roi_condition;    
catch
    roi_condition_prev = true([n_tiles, 1]); %first movement    
end
[roi_condition, logicalMask_us] = surveyedTiles(logicalMask, roiTool_metadata);

roiHandle.UserData.roi_condition = roi_condition;
roiHandle.UserData.logicalMask_us = logicalMask_us;

%Second approx.


%Loading data
if ~isequal(roi_condition, roi_condition_prev)        
    
    [labelMap_us, quantMap_cell, logMap_cell] = loadTiles( ...
    roiTool_metadata, labelMap, roi_condition);   
        
    roiHandle.UserData.labelMap_us = labelMap_us;
    roiHandle.UserData.quantMap_cell = quantMap_cell;
    roiHandle.UserData.logMap_cell = logMap_cell;
    
    % %info when crossing tile borders
    % sprintf('FoV now covers tiles: %s', ...
    %     strjoin(string(find(roi_condition)), ',')) 

else
    
    labelMap_us = roiHandle.UserData.labelMap_us;
    quantMap_cell = roiHandle.UserData.quantMap_cell;
    logMap_cell = roiHandle.UserData.logMap_cell;
end

[roiMineralMask, roiData_quant, roiData_log_n] = fuseUnfoldedTiles( ...
    logicalMask_us, labelMap_us, quantMap_cell, logMap_cell);

%% Pie chart 

pixelPopulations = zeros(1, n_masks);
for i = 1:n_masks
    pixelPopulations(i) = sum(roiMineralMask == i);
end

roiHandle.UserData.pixelPopulations = pixelPopulations;
roiHandle.UserData.roiDataXYZ = roiData_quant;
roiHandle.UserData.roiMineralMask = roiMineralMask;

%% Null matrix
%matrix: minerals x variables

temp4 = cell(1, n_masks);
for i = 1:n_masks
   temp1 = roiData_quant(roiMineralMask == i, :);
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

%% Multivariate statistical analysis 

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
       
    maskXYZ = roiData_log_n(logical, :);
    
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
