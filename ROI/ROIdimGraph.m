%% Scree and Contribution plots

function ROIdimGraph(roiHandle, choice, dataLabels, dimRange)

minerals = roiHandle.UserData.minerals; %all available
pixelNumber = roiHandle.UserData.pixelNumber; %pre-selected list
choices = roiHandle.UserData.choices;
coefs = roiHandle.UserData.coefs;
S = roiHandle.UserData.S; %eigen-values

fontSize = 8;
titleFigure = strcat('Multivariate assessment:', {' '}, minerals{choice});

idx_choice = (choices == choice);
if pixelNumber{idx_choice} == 0 
    disp('The target mineral was not intercepted.')
end

coefs = coefs{idx_choice};
S = S{idx_choice};
cos2 = coefs.^2; %matrix: variables x PC
explainedPct = 100*S./sum(S);

dimTotal = size(coefs, 2); %number of variables (= PCs)
nlab = length(dataLabels); %same as size(coefs, 2);

%understandable plot
if dimTotal>15 
    dimTotal = 15; %maximum (editable)
end
lim = find(explainedPct == 0, 1, 'first'); %only until zero %
if isempty(lim)
    lim = dimTotal; 
end

try
    contribution = 100*cos2(:, dimRange)*S(dimRange)/sum(S(dimRange)); %contribution of variable to selected PCs
catch
    contribution = zeros(size(cos2, 1), 1);
end

%% Plot 

figure(8); 

clf('reset')
hFig = gcf;
pos = get(hFig, 'Position');
set(hFig, 'Position', pos);

if dimTotal<3
    fprintf(['Insufficient variables for 3D biplots. Please, select other mineral mask or start eliminating analysis (rows) with NaN in maskXYZ \n']);
    return;
else
    
    %% Scree plot

    subplot(2, 1, 1); %Scree Plot
    
    bar(explainedPct(1:lim), 'FaceColor', [0.6463, 0.7456, 0.4188]); %light brown color
    hold on
    plot([1:lim], explainedPct(1:lim), 'k--o', 'LineWidth', 1, 'MarkerSize', 5, ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none'); 
    hold off
    %Show X and Y coordinates
    text([1:lim]', explainedPct(1:lim), strcat(num2str(explainedPct(1:lim), 2)), ...
         'horiz', 'center', 'vert', 'bottom', 'Color', 'b', 'FontSize', fontSize); 
    %without x' it was not working
    
    xlabel('PCs - eigenvalues'); 
    ylabel('%'); 
    title('Scree Plot: explained variances', 'Interpreter','none');
    
    %% Contribution plot
    
    subplot(2, 1, 2); %Contribution Plot
    
    bar(contribution, 'g'); 
    hold on
    yline(100/nlab, '--r'); %"expected average" contrib. reference line
    hold off
    
    text(nlab/2, 100/nlab+0.5, 'expected avg', 'Color', 'r', 'FontSize', fontSize); 
    
    %Show labels on top of bars. It does not work without '
    text([1:nlab]', contribution, strcat(num2str(contribution, 2)), ... 
         'horiz', 'center', 'vert', 'bottom', 'Color', 'b', 'FontSize', 0.9*fontSize); 
    
    set(gca, 'XTick', 1:nlab, 'xticklabel', dataLabels, 'fontsize', fontSize);
    xtickangle(90);
    xlabel('Variables'); 
    ylabel('%'); 
    title(strcat('Contribution to PCs: ', num2str(dimRange)));
    
    t = sgtitle(titleFigure, 'interpreter', 'none');
    t.FontSize = 1.3*fontSize; 
    t.FontWeight = 'bold';

end

end