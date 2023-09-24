close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);

DATA_DIRECTORY = fullfile(MAIN_DIRECTORY, '_data');
DTM_DIRECOTRY  = fullfile(DATA_DIRECTORY, 'ExportDMR');
PC_DIRECTORY   = fullfile(DATA_DIRECTORY, 'ExportMB');

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
view(3);
daspect([1 1 1]) % cim je posledne cislo blizsie ku 0, tym viac sa skaluje "z" suradnica
axis normal

%% DATA READING - PC
cd(MAIN_DIRECTORY);

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

%% POINTS SELECTION
ground = struct;

x = []; y = []; z = [];

for i = 1:lasCount
	x = [x; ptCloud{i}.Location(ptAttributes{i}.Classification == 2, 1)];
	y = [y; ptCloud{i}.Location(ptAttributes{i}.Classification == 2, 2)];
	z = [z; ptCloud{i}.Location(ptAttributes{i}.Classification == 2, 3)];

end

ground.Location = [x, y, z];

pcshow(ground.Location, COLOR_MED_VEG)




