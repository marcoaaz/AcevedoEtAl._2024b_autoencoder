%% Pie charts

function ROIpieChart2(roiHandle, roiTool_metadata, plotOptions)

%Input
%pixelPopulations: number of pixels per mask after registration/downscaling
%plotOptions = [1, 2] (includes both area & weight pct)

pixelPopulations = roiHandle.UserData.pixelPopulations;
density = roiTool_metadata.density; %gr/cc
minerals = roiTool_metadata.minerals; %labelled masks
triplet = roiTool_metadata.triplet; %colour, database or random

fontSize = 13;

[rankedPopulations, I] = sort(pixelPopulations, 2, 'descend');
rankedTriplet = triplet(I, :);
rankedMinerals = minerals(I);
rankedDensity = density(I);

top = rankedPopulations; %pie 1
bottom = rankedPopulations.*rankedDensity'; %pie 2

%npie: How many phases is it comfortable to watch in a pie chart? 
npie = zeros(1, 2);
top_pct = top/sum(top);
bottom_pct = bottom/sum(bottom, "omitnan");
npie(1) = sum(top_pct > 0.01); 
npie(2) = sum(bottom_pct > 0.01);

pieData = {top, bottom(~isnan(bottom))}; %inmune to NaN density row
pieTriplet = {rankedTriplet, rankedTriplet(~isnan(bottom), :)};
pieMinerals = {rankedMinerals, rankedMinerals(~isnan(bottom))};
pieTitle = {'Area %', 'Weight %'};

figure(10);
hFig = gcf;
clf('reset') %required for clearing the automatic plot

pos = get(hFig, 'Position');
set(hFig, 'Position', pos);
ha = tight_subplot(1, length(plotOptions), ...
    [.08 .1], [.03 .12], [.03 .03]);
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

for j = plotOptions
    axes(ha(j)); 

    p = pie(pieData{j});
    temp_mineral = pieMinerals{j};

    %adjusting colors of slices (issue: colormaps overwrites)
    for i = 1:length(temp_mineral)
       set(p(2*i-1), 'FaceColor', pieTriplet{j}(i, :)) 
    end
    %adjusting labels
    pText = findobj(p, 'Type', 'text'); 
    percentValues = get(pText, 'String'); 
    
    linebreak = cell(length(temp_mineral), 1);
    linebreak(:) = {newline};
    combinedtxt = strcat(temp_mineral, {':'}, linebreak, percentValues);
    for i = 1:length(temp_mineral)
        if i <= npie(j)
        pText(i).String = combinedtxt(i);
        else
        pText(i).String = "";
        end
    end

    titleH = title(pieTitle{j});    
    pos = get(titleH, 'position');
    set(titleH, 'position', [pos(1)-0.9 pos(2)-0.1])
    set(pText, 'FontSize', 0.8*fontSize, 'interpreter', 'none');
end
t = sgtitle('Abundance pie chart', 'interpreter', 'none');
t.FontSize = fontSize; 

end
