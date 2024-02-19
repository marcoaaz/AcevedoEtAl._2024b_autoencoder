clear 
clc

workingDir = 'C:\Users\n10832084\OneDrive - Queensland University of Technology\Desktop\GeoPIXE meeting';

% file1 = '65966-q1_FULL.csv';%harzburgite
% regionNames = {'clinopyroxene', 'orthopyroxene', 'olivine', 'serpentine', 'flogopite', 'fracture', 'holes', 'ilmenite', 'mineral1'}; %harzburgite

%Importing

%Input: sample
file1 = '65966-q1_FULL.csv'; %spectra (raw)
file2 = '65966-q1_imageRegionTable.csv'; %image region table (wt.%)
energyExp = 18.5; %harzburgite
desired_textX = 3.8; %check manually, default = energyExp/2
k_clusters = 8;

regionsPCT = readtable(fullfile(workingDir, file2), 'VariableNamingRule','preserve'); %values in ppm
regionNames = regionsPCT.Note; 
pixelPopulations = regionsPCT.Pixels; 
regionsGeopixe = readtable(fullfile(workingDir, file1));

%Input: Plot settings
alphaVal = 0.5;
lineWidth = 2;
fontSize = 10;
step_val = 4; %smoothing filter
percentileTH = 0.05; %plot
extra_margin = 0.05; %top and bottom

%Default: script 
n_spectra = length(regionNames);

%Referential: Line identification (GeoPIXE list)
lines_elements = {'Mg Ka', 'Al Ka', 'S Ka', 'Fe Ka', 'Mn Kb1'};
lines_energy = [1.253, 1.486, 2.307, 6.399, 6.489];
n_lines = length(lines_elements);
% 1.253, Mg, Ka_, 0.9872, K-L2, 3
% 1.486, Al, Ka_, 0.9868, K-L2, 3
% 2.307, S, Ka_, 0.9443, K-L2, 3
% 6.399, Fe, Ka_, 0.8833, K-L2, 3
% 6.489, Mn, Kb1, 0.1149, K-M3
rangeEnergy = [lines_energy(2), energyExp];

%Formatting data
x = regionsGeopixe{:, 1};
y_mtx = regionsGeopixe{:, 2:end}; 
y_mtx_normalized = y_mtx ./ sum(y_mtx, 1);
y_mtx_n_smooth = smoothdata(y_mtx_normalized, 1, 'movmean', [step_val, step_val]);
idx_rangeX = (x >= rangeEnergy(1)) & (x <= rangeEnergy(2));
y_mtx_n_smooth_crop = y_mtx_n_smooth(idx_rangeX, :);

%config axis (and text in plots)
temp_textY = max(y_mtx_n_smooth_crop, [], 'all');
temp_P = prctile(y_mtx_n_smooth_crop, [percentileTH, 100 - percentileTH], "all");

desired_P_low = 10^(log10(temp_P(1)) - extra_margin);
desired_P_high = 10^(log10(temp_P(2)) + extra_margin);
desired_textY = 10^(log10(temp_textY) - extra_margin); %check manually (text energy lines)

rangeY = [desired_P_low, desired_P_high];

[~, idx_min] = min(abs(x - desired_textX));%text location
x_min = x(idx_min);
y_min = y_mtx_n_smooth(idx_min, :); %normalizing to 1

%k-means
X = y_mtx_normalized'; %minerals x energy levels (variables)
opts = statset('Display', 'final');
[idx, C] = kmeans(X, k_clusters, 'Distance', 'cityblock',...
    'Replicates', 5, 'Options', opts);

clusteredMinerals = cell(1, k_clusters);
for j = 1:k_clusters
    temp_ind = (idx == j);
    clusteredMinerals{j} = regionNames(temp_ind);
end

% Plot 
close all

%Plot individually
cmap = [colormap(jet(n_spectra)), repmat(alphaVal, n_spectra, 1)];

hFig = figure;
hFig.Position = [100, 100, 1400, 650];

for i = 1:n_spectra        
    y = y_mtx_n_smooth(:, i);

    h(i) = plot(x, y, 'Color', cmap(i, :), 'LineWidth', lineWidth, 'Marker', 'none', 'LineStyle','-');
    hold on
end

%element vertical lines
xline(energyExp, 'Color', 'k', 'LineWidth', lineWidth)
for k = 1:n_lines
    xline(lines_energy(k), 'Color', 'k', 'LineWidth', lineWidth, 'LineStyle','--')
end
text([lines_energy, energyExp], repmat(desired_textY, 1, n_lines+1), ...
    [lines_elements, 'energyExp'], 'Interpreter','none', ...
    'FontSize', fontSize, 'HorizontalAlignment','center', 'Color', 'blue')
%mineral names
text(x_min*ones(1, n_spectra), y_min, regionNames, ...
    'Interpreter','none', 'FontSize', 0.7*fontSize, 'Color', [.5, .5, .5])
hold off

%Customize
grid on
set(gca, 'YScale', 'log')
xlim(rangeEnergy)
ylim(rangeY)

xlabel('Energy (KeV)')
title('Normalized spectra: [Al Ka, energyExp]')
legend(h, regionNames, 'Location', 'eastoutside', 'Interpreter','none', 'FontSize', fontSize)


%Plot clustered
cmap_clustered = colormap(jet(k_clusters));

hFig = figure;
hFig.Position = [100, 100, 1400, 650];

for i = 1:n_spectra    
    y = y_mtx_n_smooth(:, i);
    
    m = idx(i);

    h(i) = plot(x, y, 'Color', cmap_clustered(m, :), 'LineWidth', lineWidth, 'Marker', 'none', 'LineStyle','-');
    hold on
end
xline(energyExp, 'Color', 'k', 'LineWidth', lineWidth)
%mineral names
text(x_min*ones(1, n_spectra), y_min, regionNames, ...
    'Interpreter','none', 'FontSize', 0.7*fontSize, 'Color', [.5, .5, .5])
hold off

grid on
set(gca, 'YScale', 'log')
xlim(rangeEnergy)
ylim(rangeY)

xlabel('Energy (KeV)')
title(sprintf('Clustered Normalized spectra: %d clusters', k_clusters), 'Interpreter', 'none')
legend(h, regionNames, 'Location', 'eastoutside', 'Interpreter','none', 'FontSize', fontSize)

%%
%Kaersutite (worst spectra)
x_test = x;
y_test = y_mtx_normalized(:, 29);
plot(x_test, y_test)
figure, plot(x_test, y_test)
