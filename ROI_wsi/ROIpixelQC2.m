%% Stacked bar plot (ROI interrogation)

function ROIpixelQC2(roiHandle, roiTool_metadata)

%roiMask = roiHandle.UserData.logicalMask;
%ns = sum(roiMask);

roiData_quant = roiHandle.UserData.roiDataXYZ;
labelsModality = roiTool_metadata.variableNames2;

ns = size(roiData_quant, 1);

columnSum = nansum(roiData_quant, 1); %excluding NaN values
[~, sortOrder] = sort(columnSum, 'descend');

%ranked by descending element abundance
dataSorted = roiData_quant(:, sortOrder); 
labelsSorted = labelsModality(sortOrder);

rng(10);
myColors = hsv(length(labelsSorted)); %editing colorbar 
rColors = myColors(randperm(size(myColors, 1)), :);
 
figure(60); %stacked barplot

clf('reset')
hFig=gcf;
pos = get(hFig, 'Position');
set(hFig, 'Position', pos);
%set(hFig, 'Position', [100 100 700 400]); %left-bottom width height

p = bar(1:ns, dataSorted, 0.5, 'stack');

%Formatting
axis([0 ns+1 0 inf]);
title('QuantMap QC');
xlabel('Number of pixels');
ylabel('wt.%'); 
for j = 1:length(labelsSorted)
    %    colorSet = [colorSet myColors];
    p(j).FaceColor = 'flat';
    p(j).CData = repmat(rColors(j, :), ns, 1);
end
lgd = legend(p, labelsSorted, 'Location', 'southoutside', 'NumColumns', 6);
lgd.FontSize = 8;
lgd.Title.String = 'Characteristic X-ray (wt.%)';

clear p lgd j

end
