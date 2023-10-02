close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);

% download: https://drive.google.com/drive/folders/1EWlkqu3qofB5kB5QLcPkrFyJxFm2zBec?usp=sharing
DATA_DIRECTORY = fullfile(MAIN_DIRECTORY, 'data');
DTM_DIRECOTRY  = fullfile(DATA_DIRECTORY, 'ExportDMR4');
PC_DIRECTORY   = fullfile(DATA_DIRECTORY, 'ExportMB4');

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

%%
figure("Name","Original Point Cloud")
handle.plotPtCloud3D(colorMap, [2, 3, 4, 5],"useData","original")

figure("Name","Normalized Point Cloud")
handle.plotPtCloud3D(colorMap, [2, 3, 4, 5],"useData","normalized")

%%
handle = handle.computeMetricRasters();

%%
figure
handle.plotStatRaster("plotData","Hp95")
colormap jet
colorbar
axis equal
axis xy

%% EXPORT TO TIFF
geotiffwrite("_images/data4/Hmax_data4", handle.statRasters.Hmax, rasterReference,"CoordRefSysCode",8353);

geotiffwrite("_images/data4/Hmean_data4", handle.statRasters.Hmean, rasterReference,"CoordRefSysCode",8353);

geotiffwrite("_images/data4/Hmedian_data4", handle.statRasters.Hmedian, rasterReference,"CoordRefSysCode",8353);

geotiffwrite("_images/data4/Hp25_data4", handle.statRasters.Hp25, rasterReference,"CoordRefSysCode",8353);

geotiffwrite("_images/data4/Hp75_data4", handle.statRasters.Hp75, rasterReference,"CoordRefSysCode",8353);

geotiffwrite("_images/data4/Hp95_data4", handle.statRasters.Hp95, rasterReference,"CoordRefSysCode",8353);


%%
[H, R] = readgeoraster('Hmax.tif', 'OutputType', 'double');

H = flipud(H);

