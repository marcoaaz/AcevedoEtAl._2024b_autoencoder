%% 3D Bar graph of non-zero pixels 

function ROIbarGraph3D(roiHandle, threshold, faceAlpha)

nullPixels = roiHandle.UserData.nullPixels;
fontSize = 9;
opacityVal = 1;


matrix = nullPixels(:, 2:end);
variables = matrix.Properties.VariableNames;
n_minerals = size(nullPixels, 1); %mineral axis ticks
n_variables = length(variables); %element axis ticks

figure(30);
hFig = gcf;
pos = get(hFig, 'Position');
set(hFig, 'Position', pos);
%set(hFig, 'Position', [250 100 850 600]);

h = bar3(100-matrix{:, :}); %getting handle
hold on

colormap(flipud(gray(10))) %colormap for 'interp'
set(h,'FaceAlpha', opacityVal); %Option: opacity <1 for 3D visualization
for k = 1:length(h)
    zdata = h(k).ZData;
    h(k).CData = zdata;
    h(k).FaceColor = 'interp'; %interpolated color
end
colorbar
view([0 90])
ax = gca;
ax.Interactions = [rotateInteraction dataTipInteraction];

%x-axis
xhandle = xlabel('Element');
set(gca, 'XTickLabel', variables', 'fontsize', fontSize*0.9) % Change tick labels
xticks(1:n_variables); 
xtickangle(40)
set(xhandle, 'FontSize', fontSize) 
set(xhandle,'FontWeight','bold') %bold font
xlim([0, n_variables+1]);

%y-axis
yhandle = ylabel('Mineral');
set(gca, 'YTickLabel', nullPixels.minerals, 'fontsize', fontSize*0.9)
yticks(1:n_minerals);
ytickangle(40)
set(yhandle, 'FontSize', fontSize) 
set(yhandle,'FontWeight','bold') %bold font
ylim([0, n_minerals+1]);

%z-axis
zhandle = zlabel('Percentage %');
set(zhandle, 'FontSize', fontSize) 
set(zhandle,'FontWeight','bold') %bold font

%appereance
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'none'; 
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'none'; 

if faceAlpha == 0
    title('Positive Pixels', 'FontSize', fontSize*1.2)
elseif faceAlpha > 0
    title(sprintf('Positive Pixels (threshold plane: %d pct)', threshold), 'FontSize', fontSize*1.2)
end

%Plotting threshold plane (can be LOD)
X1 = [0.6 n_variables + 0.4; 
    0.6 n_variables + 0.4];
Y1 = [0 0; 
    n_minerals+1 n_minerals+1];
Z1 = [100-threshold 100-threshold; 
    100-threshold 100-threshold];
surf(X1, Y1, Z1, 'FaceColor', 'blue', 'edgecolor', 'none', 'FaceAlpha', faceAlpha); 
%Change 'FaceAlpha' to 1 (below) to see selected variables for PCA
hold off



clear maxis eaxis h k 

end