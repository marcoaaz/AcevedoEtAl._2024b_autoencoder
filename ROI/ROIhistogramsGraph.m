%% Histograms of all channels

function ROIhistogramsGraph(roiHandle, choices, dataLabels, pctOut, xAxisLabel, mode)
%mode: "logarithmic", "linear"


roiDataXYZ = roiHandle.UserData.roiDataXYZ;
roiMineralMask = roiHandle.UserData.roiMineralMask;
minerals = roiHandle.UserData.minerals;

n_bins = 40;
fontSize = 12;
faceColor = [0.9290 0.6940 0.1250];

%Settings (big plot)
nB = length(dataLabels); %number of plots
nf = ceil(nB^0.5); %distribution of sub-plots
if nf^2 - nB >= nf
    nrows = nf-1;
    ncolumns = nf;
else 
    nrows = nf;
    ncolumns = nf;
end
clear nf

k = 0;
for j = choices
    k = k + 1;

    choice = j;

    titleFigure = strcat(' Histograms ', ' (', mode, ' scale):', {' '}, minerals(choice, 1));
    maskXYZ = roiDataXYZ(roiMineralMask == choice, :);
    
    % Log-normal/linear histograms
    figure(40 + k);
    
    clf('reset') %required for clearing the automatic plot
    hFig= gcf;
    
    pos = get(hFig, 'Position');
    set(hFig, 'Position', pos);
    %set(hFig, 'Position', [50 30 1100 700]);
    
    ha = tight_subplot(nrows, ncolumns, [.11 .07], [.1 .1], [.09 .09]);
    %[ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    
    for i = 1:nB
        axes(ha(i));
        
        temp1 = maskXYZ(:, i);
        %temp1(isnan(temp1)) = []; %exclude them, since they are in CPS not in ppm
        
        [temp1, ~, ~] = rmoutliers(temp1, 'percentiles', [0, 100-pctOut]);%Optional
    
        if mode == "linear"
            temp1 = temp1(~(temp1 == 0)); %useful for EDX maps
            %for LA maps set Iolite DRS to:
            %Threshold to use when masking low signals = 0
            %Seconds to trim before/after low signals = 0   
    
            histogram(temp1, n_bins, 'FaceColor', faceColor); 
    
        elseif mode == "logarithmic"
    
            histogram(log10(temp1), n_bins, 'FaceColor', faceColor); 
        end    
        grid on
        ax = gca;
        ax.YAxis.Exponent = 2;
        ax.XAxis.Exponent = 4; %for Synchrotron XFM =4, depends on the technique
    
        xlabel(xAxisLabel, 'FontSize', fontSize*0.8);
        %xlim([-1 inf]); 
        ylim([0 inf]); %clipped input data
        
        title(dataLabels{i});    
    end
    t = sgtitle(titleFigure, 'interpreter', 'none');
    t.FontSize = fontSize; 
end

end