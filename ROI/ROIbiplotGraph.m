%% Tridimensional biplot of eigen vectors

function ROIbiplotGraph(roiHandle, choices, dataLabels, dimRange)

coefs = roiHandle.UserData.coefs;
S = roiHandle.UserData.S;
maskLabel = roiHandle.UserData.maskLabel; %minerals(choices, 1);

n_choices = length(choices);
fontSize = 10;

titleFigure = strcat('Composition eigenvectors in PC: ', num2str(dimRange));

%Settings (big plot)
nB = n_choices; %number of plots
nf = ceil(nB^0.5); %distribution of sub-plots
if nf^2 - nB >= nf
    nrows = nf-1;
    ncolumns = nf;
else 
    nrows = nf;
    ncolumns = nf;
end
clear nf

figure(67);

clf('reset') %required for clearing the automatic plot
hFig=gcf;

pos = get(hFig, 'Position');
set(hFig, 'Position', pos);
%set(hFig, 'Position', [50 30 1100 700]);
ha = tight_subplot(nrows, ncolumns, [.12 .1], [.08 .1], [.1 .08]);
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

for i = 1:n_choices
    
    coefs_temp = coefs{i};

    obsTotal = size(coefs_temp, 1); %number of eigen-vectors
    dimTotal = size(coefs_temp, 2); %number of variables (= PCs)
    cos2 = coefs_temp.^2; 
    
    if dimTotal > 0
        quality = 100*sum(cos2(:, dimRange), 2); %quality of representation in factor map (biplot)
        explainedPct = 100*S{i}./sum(S{i}); %pct of explained variance 
    
        %Colormap for biplot
        map = zeros(obsTotal, 3); %each row being a RGB triplet
                
        map(:, 1) = quality./max(quality); %R
        map(:, 2) = 0; %G
        map(:, 3) = 1-quality./max(quality); %B
        
        %figure
        axes(ha(i));
    
        h = biplot(coefs_temp(:, dimRange), 'Color', 'b', 'VarLabels', dataLabels); 
        
        %Layout
        xlabel(strcat('PC', num2str(dimRange(1)), sprintf('= %0.1f pct', ...
            round(explainedPct(dimRange(1)), 2))));
        ylabel(strcat('PC', num2str(dimRange(2)), sprintf('= %0.1f pct', ...
            round(explainedPct(dimRange(2)), 2))));
        if length(dimRange)>2
            zlabel(strcat('PC', num2str(dimRange(3)), sprintf('= %0.1f pct', ...
            round(explainedPct(dimRange(3)), 2))));
        end
        title(maskLabel{i}, 'FontWeight', 'normal', ...
            'FontSize', fontSize, 'interpreter', 'none');
    
        %colouring eigenvectors
        %lines
        for k = 1:obsTotal
            h(k).Color = map(k, :); 
            h(k).LineStyle = '-'; 
        end
        %obs.
        for k = (obsTotal+1):(2*obsTotal)
            h(k).MarkerEdgeColor = 'k';  %black
            h(k).MarkerSize = 10; 
        end
        %text
        for k = (2*obsTotal+1):(3*obsTotal)
            h(k).FontSize = 0.7*fontSize;  % Specify character size
            h(k).FontWeight = 'bold';  % Specify bold font
            h(k).Color = 'k';  
        end

    else
        continue
    end
end
t = sgtitle(titleFigure);
t.FontSize = 1.1*fontSize; 
t.FontWeight = 'bold';

end