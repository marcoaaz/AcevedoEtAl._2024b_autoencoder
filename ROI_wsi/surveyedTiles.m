function [roi_condition, logicalMask_us] = surveyedTiles( ...
    logicalMask, roiTool_metadata)

%Input
%localMask = ROI as a binary mask with the same size as 'ref_img'. Note: This
%could be improved to read polygon types vertex and save RAM memory
%scale_factor = upscale localMask to the size of the stack (much larger)
%corners_array = top-left and bottom-right pixel coordinates of each tile

%Output
%roi_condition = array with the image tiles of interest (ToI)

scale_factor = roiTool_metadata.scale_factor;
corners_array_WH = roiTool_metadata.corners_array_WH;
corners_array_WH_ds = roiTool_metadata.corners_array_WH_ds;

%% Tiles containing FoV

s = regionprops(logicalMask, 'BoundingBox');
bbox = s.BoundingBox;
bbox2 = bbox*scale_factor; %upscaling

bbox_tl_col = bbox2(1); %x
bbox_tl_row = bbox2(2); %y
bbox_br_col = bbox_tl_col + bbox2(3); %W
bbox_br_row = bbox_tl_row + bbox2(4); %H

%first approx.
tl_condition = ((corners_array_WH(:, 3) > bbox_tl_row) & ...
    (corners_array_WH(:, 4) > bbox_tl_col));
br_condition = ((corners_array_WH(:, 1) < bbox_br_row) & ...
    (corners_array_WH(:, 2) < bbox_br_col));
roi_condition = tl_condition & br_condition; %index vector
n_load = sum(roi_condition, 'all');

%% Surveyed region of interest (ROI): upscaled and tiled

tile_ind = 1:size(corners_array_WH, 1);
tile_ind2 = tile_ind(roi_condition);
corners_array2 = corners_array_WH(roi_condition, :);
corners_array2_ds = corners_array_WH_ds(roi_condition, :); %downscaled sub

logicalMask_us = cell(1, n_load); %upscaled tiles

parfor k = 1:n_load %depends on img_size          
    %sprintf('tile %.0f has %.0f rows and %.0f cols', tile_ind2(k), n_rows, n_cols)
    
    %extracting tiles
    temp_logical = logicalMask( ...
        corners_array2_ds(k, 1):corners_array2_ds(k, 3), ...
        corners_array2_ds(k, 2):corners_array2_ds(k, 4));
    
    temp_logical2 = imresize(temp_logical, ...
        [corners_array2(k, 5), corners_array2(k, 6)], "nearest");
               
    %storing in cache
    logicalMask_us{k} = temp_logical2;    
end

end