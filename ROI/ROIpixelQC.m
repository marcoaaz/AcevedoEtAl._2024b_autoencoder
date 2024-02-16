%% Stacked bar plot (ROI interrogation)

function ROIpixelQC(roiHandle, xyzModality, labelsModality)

% probe (ROI)
roiMask = roiHandle.UserData.logicalMask;
ns = sum(roiMask);

columnSum = nansum(xyzModality, 1); %excluding NaN values
[~, sortOrder] = sort(columnSum, 'descend');
dataSorted = xyzModality(:, sortOrder); %ranked by descending element abundance
labelsSorted = labelsModality(1, sortOrder);

rng(10);
myColors = hsv(length(labelsSorted)); %editing colorbar 
rColors = myColors(randperm(size(myColors, 1)), :);
 
figure(60); %stacked barplot

clf('reset')
hFig=gcf;
pos = get(hFig, 'Position');
set(hFig, 'Position', pos);
%set(hFig, 'Position', [100 100 700 400]); %left-bottom width height

C = dataSorted(roiMask, :);
p = bar(1:ns, C, 0.5, 'stack');
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
