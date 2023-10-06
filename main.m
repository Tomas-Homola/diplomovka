close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);

% download: https://drive.google.com/drive/folders/1EWlkqu3qofB5kB5QLcPkrFyJxFm2zBec?usp=sharing
DATA_DIRECTORY = fullfile(MAIN_DIRECTORY, 'data');
DTM_DIRECOTRY  = fullfile(DATA_DIRECTORY, 'ExportDMR_lesBodiky');
PC_DIRECTORY   = fullfile(DATA_DIRECTORY, 'ExportMB_lesBodiky');

cd(MAIN_DIRECTORY);

%% COLOR MAP
colorMap = load('las_colormap.txt')/255; % color map for classified point cloud

COLOR_GROUND     = colorMap(2,:);
COLOR_VEG_LOW    = colorMap(3,:);
COLOR_VEG_MEDIUM = colorMap(4,:);
COLOR_VEG_HIGH   = colorMap(5,:);

%% DATA READING - DMR
cd(DTM_DIRECOTRY);

[heightMatrix, rasterReference] = readgeoraster('dmr.tif', 'OutputType', 'double');
infoDMR = georasterinfo('dmr.tif');

heightMatrix = standardizeMissing(heightMatrix, infoDMR.MissingDataIndicator);
heightMatrix = flipud(heightMatrix);
rasterReference.ColumnsStartFrom = 'south';

%% PLOT - DMR
figure('Name', 'Digital Terrain Model data')
title 'DTM data'
mapshow(heightMatrix, rasterReference, 'DisplayType', 'surface');
demcmap(heightMatrix, 256)
view(2);
daspect([1 1 1]) % cim je posledne cislo blizsie ku 0, tym viac sa skaluje "z" suradnica
axis normal

%% DATA READING - PC
cd(PC_DIRECTORY);
lasFiles = dir('*.las'); % find all *.las files
lasCount = length(lasFiles); % number of found *.las files

ptCloud = cell(lasCount, 1);
ptAttributes = cell(lasCount, 1);
colorData = cell(lasCount, 1);

for i = 1:lasCount
    lasReader = lasFileReader(lasFiles(i).name);
    [ptCloud{i}, ptAttributes{i}] = readPointCloud(lasReader, 'Attributes', 'Classification');
    colorData{i} = reshape(label2rgb(ptAttributes{i}.Classification, colorMap, 'k'), [], 3);
end

%% PLOT - PC
% 2 -> Ground
% 3 -> Low Vegetation
% 4 -> Medium Vegetation
% 5 -> High Vegetation
selectedClasses = [2, 3, 4, 5];
	
figure('Name', 'Point Cloud data')
title 'Point Cloud data'
hold on
for i = 1:lasCount
	classMember = ismember(ptAttributes{i}.Classification, selectedClasses);
% 	classMember = classMember & ptCloud{i}.Intensity < 10;
		if any(classMember)
			pcshow(ptCloud{i}.Location(classMember, :), colorData{i}(classMember, :), ...
				'MarkerSize', 20)
		end
end
hold off

%%
cd(MAIN_DIRECTORY);

handle = dataHandler(heightMatrix, rasterReference, ptCloud, ptAttributes);
handle = handle.meshPlane();
handle = handle.computePointCloudAttributes();
handle = handle.normalizePtCloud("method","DTM");

%% PLOT ORIGINAL AND NORMALIZED POINT CLOUD
figure("Name","Original Point Cloud")
handle.plotPtCloud3D(colorMap, [2, 3, 4, 5],"useData","original")

figure("Name","Normalized Point Cloud")
handle.plotPtCloud3D(colorMap, [2, 3, 4, 5],"useData","normalized")

%% COMPUTE ALL METRICS FOR GIVEN POINT CLOUD
handle = handle.computeMetricRasters();

%% PLOT SELECTED METRIC
figure
handle.plotMetricRaster("plotData","PPR")
% colormap jet
colormap hot
colorbar
axis equal
axis xy

%%
figure
imagesc([handle.x1 handle.x2], [handle.y1 handle.y2], temp,...
	'AlphaData', handle.alphaData_DTM)
title 'Max of normalized height'
colormap hot
colorbar
axis equal
axis xy

%% EXPORT SELECTED METRIC TO .TIF FILE
handle.exportMetricRaster("data_lesSipka","exportLayer","Hstd");

%% EXPORT ALL METRICS TO .TIF FILES
handle.exportAllMetricRasters("data_lesBodiky");

%%
[H, R] = readgeoraster('Hmax.tif', 'OutputType', 'double');

H = flipud(H);

%% Export color palette
c = jet;
c = c*255;
writematrix(c, "jet_palette.txt","Delimiter",' ');




