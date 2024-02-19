%% Figure # script
clear 
clc

%GeoPIXE convention
%spectra plots: green (data), red (fit), purple (bg)
%spectra exports columns: raw data of named region, SNIP= fitted background, model fit

%Chemical Element Info
% 1.253, Mg, Ka_, 0.9872, K-L2, 3
% 1.486, Al, Ka_, 0.9868, K-L2, 3
% 2.307, S, Ka_, 0.9443, K-L2, 3
% 6.399, Fe, Ka_, 0.8833, K-L2, 3
% 6.489, Mn, Kb1, 0.1149, K-M3

lines_elements = {'Mg Ka', 'Al Ka', 'S Ka', 'Fe Ka', 'Mn Kb1'};
lines_energy = [1.253, 1.486, 2.307, 6.399, 6.489];
n_lines = length(lines_elements);

utilitaryDir = 'E:\Alienware_March 22\02-Research geology\05-QUT\Paper 3\chris ryan collaboration\GeoPIXE meeting';
table_intensity = readtable(fullfile(utilitaryDir, 'energies_intensity.xlsx'), 'VariableNamingRule','preserve');

%Peak database
emissionElements2 = {
    'Th', 'Nd', ... %phosphate
    'Zr', 'Ti',... %residual    
    };
emissionTransition2 = {
    'L', 'L', ...
    'K', 'K', ...    
    };
% emissionElements2 = {
%     'Th', 'Nd', ... %phosphate
%     'Zr', 'Ti',... %residual
%     'La' %phosphate
%     };
% emissionTransition2 = {
%     'L', 'L', ...
%     'K', 'K', ...
%     'L'
%     };

n_emission2 = length(emissionElements2);
emissionColours1 = lines(n_emission2); %colormap
idx_per = randperm(n_emission2);
emissionColours2 = emissionColours1(idx_per, :);

cell_intensity = cell(1, n_emission2);
for m = 1:n_emission2
    idx_1 = strcmp(emissionElements2{m}, table_intensity.Element);
    idx_2 = contains(table_intensity.Line, emissionTransition2{m});
    idx_3 = idx_1 & idx_2;

    cell_intensity{m} = table_intensity(idx_3, :);
end


%%
workingDir = fullfile(utilitaryDir, '\MarcoSpectraE3_66039\MarcoSpectraE3_66039');

file1 = '66039_SumCal.csv'; %sum
file2 = '66039_Back_2488Pixels.csv'; %mask_back
file3 = '66039_Apatite_14386Pixels.csv'; %mask

n_pixels = [22625828, 2488, 14386]; %known from GeoPIXE export (Michael)

S = sprintf('%0.2e %0.2e %0.2e', n_pixels(1), n_pixels(2), n_pixels(3));
S = regexprep(S, {'e[+-]0+\>', 'e\+?(-?)0*(\d+)'}, {'', '{\\times}10^{$1$2}'});
S_split = strsplit(S, ' ');

%%
%Import
cd(workingDir)
table1 = readtable(file1, 'VariableNamingRule','preserve');
table2 = readtable(file2, 'VariableNamingRule','preserve');
table3 = readtable(file3, 'VariableNamingRule','preserve');

%Circumstantial info
legendName = 'Calibrated spectrum:';
cell_regionNames = {'Overall map area:', 'Residual spectra:', 'Phosphate mask:'};
for k = 1:length(cell_regionNames)
    cell_regionNames{k} = strcat(cell_regionNames{k}, {' '}, S_split(k), ' pixels');
end
cell_regionNames = string(cell_regionNames);
cell_lineStyle = {'-', ':', '-'}; %raw, background, model
cell_table = {table1, table2, table3};

data_cmap = [
    [0, 1, 1, 0.7]; 
    [1, 0, 1, 0.7]; 
    [1, 165/255, 0, 0.7]
    ]; %C, M, Y

energyExp = 18.5; %keV
rangeEnergy = [lines_energy(2), energyExp + .3];
desired_textY = 10^8;
lineWidth = 4;
lineWidths = [lineWidth/2, lineWidth, lineWidth];
fontSize = 10;
% rangeY = [0.0001, Inf];
greyValue_raw = 0.6;
minVal = 0.0001; %base of X-ray peak (database)
scaleVal = 10^5;
%% Plot 
close all

hFig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
 
for i = 1:3 %mask
    temp_table = cell_table{i};

    for j = [1, 3] %curve: raw, background, model
        x = temp_table{:, 1};
        y = temp_table{:, 1 + j};
                       
        %Data for customisation
        if j == 1 %raw
            temp_cmap = greyValue_raw*[1, 1, 1]; %grey

            %for normalisation
            [~, temp_idx] = min(abs(x - energyExp));
            y_val = y(temp_idx);
        else
            temp_cmap = data_cmap(i, :);
        end

        %Normalisation (to elastic scatter at 18.5 keV)
        y_norm = y/y_val;
        % y_norm = y;

        %label text (elements)
        if (i == 1) && (j == 1) %manually selected: all, raw
            x_temp = x;
            y_temp = y_norm;            
        end        
        
        %Data for plot
        h(i) = plot(x, y_norm, 'LineWidth', lineWidths(j), 'Marker', 'none', ...
            'LineStyle', cell_lineStyle{j}, 'Color', temp_cmap);        
        hold on
        
    end    
end

%Plot customisation

%element vertical lines
emissionEnergies2 = zeros(1, n_emission2);
correspondingY2 = zeros(1, n_emission2);
correspondingY2_raw = zeros(1, n_emission2);
for s = 1:n_emission2
    intensityTable_temp = cell_intensity{s};
    energy_temp = (intensityTable_temp.Energy_eV)/1000;
    intensity_temp = (intensityTable_temp.('Relative intensity'))/scaleVal;
    
    lineColour = emissionColours2(s, :);
    n_peaks = length(intensity_temp);
    for w= 1:n_peaks
        energy_temp2 = energy_temp(w);               

        h_lines(w) = plot(energy_temp2*[1, 1], [minVal, minVal + intensity_temp(w)], ...
            'LineWidth', lineWidth/2, 'Color', lineColour, 'LineStyle', '-');

        %y-coordinates (element labels)
        % emissionEnergies2(s) = mean(energy_temp2, "all"); %option 1: avg. emission energies                
    end
    %option 2: max peak emission energies (on top of raw curve) 
    [intensity_max, idx_max] = max(intensity_temp, [], "all");       
    emissionEnergies2(s) = energy_temp(idx_max);
    [~, temp_idx] = min(abs(x_temp - emissionEnergies2(s))); %text Y value
    correspondingY2_raw(s) = y_temp(temp_idx);
    
    %option 3: peak intensity value    
    correspondingY2(s) = intensity_max; 
end
text(emissionEnergies2, correspondingY2, ...
    emissionElements2, 'Interpreter','none', 'FontSize', fontSize*1.3, ...
    'HorizontalAlignment','center', 'Color', 'black', ...
    'Rotation', 45, 'FontWeight','bold');
%annotation('textarrow',x,y,'String',' Growth ','FontSize',13,'Linewidth',2)
hold off
axis on

%Customize
xlim(rangeEnergy)
% ylim(rangeY)
grid on
grid(gca,'minor')

ax = gca;
set(ax, 'YScale', 'log')
set(ax, 'YMinorGrid','off')
ax.XAxis.FontSize = fontSize*1.5;
ax.YAxis.FontSize = fontSize*1.5;


xlabel('Energy (KeV)', 'FontSize', fontSize*2)
ylabel('Normalised counts', 'FontSize', fontSize*2)
title('GeoPIXE spectra', 'FontSize', fontSize*2)
leg = legend(h, cell_regionNames, 'Location', 'southoutside', ...
    'FontSize', fontSize*1.5, 'Orientation', 'horizontal');
title(leg, legendName, 'FontSize', fontSize*1.5)

