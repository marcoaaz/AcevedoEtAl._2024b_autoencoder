
function [labelMap_us, quantMap_cell, logMap_cell] = loadTiles( ...
    roiTool_metadata, labelMap, roi_condition)

%Input
%original_path = path of image tile stack pyramid containing original pixels
%tileNames_sorted = name of tiles in original_path
%labelMap = segmentation map from QuPath
%logicalMask = ROI Tool mask in low resolution ('ref_img')
%scale_factor = multiply for upscaling
%good_idx = chemical elements to choose from

%Output
%tiles_unfolded = appended data for ROITool interrogation

original_path = roiTool_metadata.original_path;
tileNames_sorted= roiTool_metadata.tileNames_sorted;
corners_array_WH = roiTool_metadata.corners_array_WH;
corners_array_WH_ds = roiTool_metadata.corners_array_WH_ds;
good_idx = roiTool_metadata.good_idx;


%%

n_load = sum(roi_condition, 'all');
tileNames_sorted2 = tileNames_sorted(roi_condition);

corners_array2 = corners_array_WH(roi_condition, :);
corners_array2_ds = corners_array_WH_ds(roi_condition, :); %downscaled sub

%Loading ROI mask and segmentation (tiles)
labelMap_us = cell(1, n_load); %upscaled tiles
filePaths = cell(1, n_load);

parfor k = 1:n_load %depends on img_size  
   
    %extracting tiles
    temp_label = labelMap( ...
        corners_array2_ds(k, 1):corners_array2_ds(k, 3), ...
        corners_array2_ds(k, 2):corners_array2_ds(k, 4));
    
    temp_label2 = imresize(temp_label, ...
        [corners_array2(k, 5), corners_array2(k, 6)], "nearest");
           
    %storing in cache    
    labelMap_us{k} = temp_label2;
    filePaths{k} = fullfile(original_path, tileNames_sorted2{k});%paths
end

%Loading quantitative image tile stack

quantMap_cell = cell(1, n_load); %upscaled tiles
logMap_cell = cell(1, n_load); %upscaled tiles

parfor m = 1:n_load

    temp_tile = imread(filePaths{m});
    
    %medicine: clip zeros (for meaningful histograms and plots)
    quantMap = temp_tile;
    quantMap(temp_tile < 0) = 0;
    
    %Force positive distribution (including negatives reduces noise)
    temp1 = temp_tile - min(temp_tile, [], 'all'); 
    
    %log-transformation and preventing NaN (for meaningful dendrograms)
    logMap = log(temp1 + 1);
       
    %storing in cache
    quantMap_cell{m} = quantMap(:, :, good_idx);
    logMap_cell{m} = logMap(:, :, good_idx);
end

end

