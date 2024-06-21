%% Dendrogram from Hierarchical Cluster Analysis of Principal Components Analysis eigenvectors

function ROIdendrogramPCA2(roiHandle, roiTool_metadata, ...
    dimRange, inner, distance)

%inner = 'ward'; %inner squared distance (method)

maskLabel = roiHandle.UserData.maskLabel;
pixelNumber = roiHandle.UserData.pixelNumber;
coefs = roiHandle.UserData.coefs;
choices = roiTool_metadata.choices;
dataLabels = roiTool_metadata.variableNames2;

fontSize = 10;

%Settings (small plot)
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
figure(49);

clf('reset') %required for clearing the automatic plot
hFig = gcf;

pos = get(hFig, 'Position');
set(hFig, 'Position', pos);

%set(hFig, 'Position', [50 30 1100 700]);
ha = tight_subplot(nrows, ncolumns, [.14 .12], [.12 .12], [.1 .1]);
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%vertical, horizontal

for i = 1:n_choices
    
    coefs_temp = coefs{i};
    dimTotal = size(coefs_temp, 2);    
    if dimTotal >= dimRange %0
        Ypd = pdist(coefs_temp(:, dimRange)); %euclidean distance
        %Ysf = squareform(Ypd); %(nxn)
    
        Zln = linkage(Ypd, inner); %inner = Ward linkage (cluster squared distances)
        leafOrder = optimalleaforder(Zln, Ypd);
        
        % topVal = max(Zln(:, 3))/3; %x-axis plotting range
        topVal = Inf;

        %correlation coef. between distances (~1: clustering reflects your data well)
        co = cophenet(Zln, Ypd);         
               
        axes(ha(i));
        
        %'Criterion', 'inconsistency' --> doesn't work well for coloring dendrograms
        H = dendrogram(Zln, 0, 'reorder', leafOrder, 'Orientation', 'left' ...
            , 'ColorThreshold', distance, 'Labels', dataLabels);
        xline(distance, '--r');
               
        %Layout        
        xlim([0 topVal]);
        set(H, 'LineWidth', 2);
        xlabel('Height', 'FontSize', fontSize);                      
        % ylabel('Ordered leafs (optimal)', 'FontSize', fontSize);
        
        tsmall = title( ...
            string(strcat(maskLabel{i}, {', '}, CommaFormat(pixelNumber{i}),...
            ' px, coph.=', {' '}, num2str(round(co, 2)))), ...
            'interpreter', 'none');
        tsmall.FontSize = fontSize;

        % %Optional calculation:
        % %matrix(:,4) of inconsistency coeff. 
        % I = inconsistent(Zln); 
        % %vector of respective clusters for original variables
        % %if we require to obtain dendrogram leafs given a threshold
        % T = cluster(Zln, 'cutoff', distance, 'Criterion', 'distance');
    else
        continue
    end
end
t = sgtitle(strcat('HCA Dendrogram: PC', {' '}, num2str(dimRange)));
t.FontSize = 1.2*fontSize;

end