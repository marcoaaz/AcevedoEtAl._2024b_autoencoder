%% Bivariate binscatter plots (all assemblage)

function ROIbinscatterGraph2(roiHandle, roiTool_metadata, x1, x2, pctOut)

%Input
%x1 and x2: chemical elements to study
%pctOut: percentile to crop element histograms and define general axis limits

roiDataXYZ = roiHandle.UserData.roiDataXYZ;
roiMineralMask = roiHandle.UserData.roiMineralMask;
pixelPopulations = roiHandle.UserData.pixelPopulations;
dataLabels = roiTool_metadata.variableNames2; %str of chemical elements
minerals = roiTool_metadata.minerals;
choices = roiTool_metadata.choices;

fontSize = 12;
transparencyVal = 0.8;
n_pixels = sum(pixelPopulations); %within ROI
dn1 = strcmp(dataLabels, x1); %look up for specific columns
dn2 = strcmp(dataLabels, x2); 

bottom_val = 0; %axis limit (default)
top_val = Inf;

%Generalising axis limits (using percentiles)
x_list = [];
y_list = [];
for j = choices    
    x = roiDataXYZ(roiMineralMask == j, dn1); 
    y = roiDataXYZ(roiMineralMask == j, dn2); 

    index = ~(x == 0) & ~(y == 0); %avoiding pixels at origin       
    xValue = x(index);
    yValue = y(index);   
      
    x_list = [x_list; xValue];    
    y_list = [y_list; yValue];        
end
P_x = prctile(x_list, [pctOut, 100-pctOut], "all");
P_y = prctile(y_list, [pctOut, 100-pctOut], "all");

x_bottom = max(bottom_val, P_x(1));
x_top = min(P_x(2), top_val);
y_bottom = max(bottom_val, P_y(1));
y_top = min(P_y(2), top_val);

%medicine
if x_bottom == x_top
    x_top = x_top + 1;
end
if y_bottom == y_top
    y_top = y_top + 1;
end

%Plot grid settings (small plot)
n_choices = length(choices);
nB = n_choices; %number of plots
nf = ceil(nB^0.5); %distribution of sub-plots
if nf^2 - nB >= nf
    nrows = nf-1;
    ncolumns = nf;
else 
    nrows = nf;
    ncolumns = nf;
end
clear nf nB

%Plot
figure(50);

clf('reset') %required for clearing the automatic plot
hFig=gcf;

pos = get(hFig, 'Position');
set(hFig, 'Position', pos);
%set(hFig, 'Position', [50 30 1100 700]);

ha = tight_subplot(nrows, ncolumns, [.14 .07], [.14 .16], [.12 .06]);
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w) 
%vertical, horizontal

k = 0;
for i = choices
    k = k+1;   
    axes(ha(k));
    
    %filters
    x = roiDataXYZ(roiMineralMask == i, dn1); 
    y = roiDataXYZ(roiMineralMask == i, dn2); 

    index = ~(x == 0) & ~(y == 0); %avoiding pixels at origin       
    xValue = x(index);
    yValue = y(index);       

    %plot
    if ~isempty(xValue) && ~isempty(yValue)
       L(i) = binscatter(xValue, yValue);
    else
       L(i) = binscatter([1 2], [1 2]);
    end
    L(i).FaceAlpha = transparencyVal;
    
    grid on
    ax = gca;
    
    %axis
    ax.YAxis.Exponent = 4;
    ax.XAxis.Exponent = 4; %for Synchrotron XFM =4, depends on the technique

    set(gca,'Layer','top','GridColor','k','GridAlpha', 1)
    colormap(gca, 'parula'); 
    colorbar('FontSize', 0.6*fontSize);    
    
    xlabel(x1); 
    ylabel(x2); %strcat(x2, {' '}, '%')
    xlim([x_bottom x_top]);  %xtop
    ylim([y_bottom y_top]); %ytop

    % xticks(0:10:xtop);
    % yticks(0:10:ytop);    
    % xa = get(gca, 'XTickLabel');  
    % yb = get(gca,'YTickLabel');  
    % set(gca, 'XTickLabel', xa, 'fontsize', 7)      
    % set(gca,'YTickLabel', yb,'fontsize',7)  
    
    title([sprintf(strcat(minerals{i}, ' = %s px'), CommaFormat(pixelPopulations(i))), newline], ...
        'interpreter', 'none');   
end
hold off
t = sgtitle(sprintf('Binscatter plots (ROI = %.f px)', n_pixels));
t.FontSize = fontSize;

end