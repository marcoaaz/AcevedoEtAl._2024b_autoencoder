function [roiMineralMask, roiData_quant, roiData_log_n] = fuseUnfoldedTiles( ...
    logicalMask_us, labelMap_us, quantMap_cell, logMap_cell)

n_load = length(logicalMask_us);

%Sequential unfolding (matrix --> vector; XYZ generator)
roiMineralMask = [];
roiData_quant = []; %for geochemical quantification
roiData_log = []; %for multi-variate statistics
for n = 1:n_load
    temp_logical = logicalMask_us{n}(:);
    n_pixels_temp = length(temp_logical);
    
    temp_idx = reshape(labelMap_us{n}, n_pixels_temp, 1);%mineral indexes
    temp_quant = reshape(quantMap_cell{n}, n_pixels_temp, []);
    temp_log = reshape(logMap_cell{n}, n_pixels_temp, []);
    
    %appending
    roiMineralMask = [roiMineralMask; temp_idx(temp_logical)];
    roiData_quant = [roiData_quant; temp_quant(temp_logical, :)];
    roiData_log = [roiData_log; temp_log(temp_logical, :)];

end

roiData_log_n = normalize(roiData_log, 1, 'range');

end