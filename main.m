close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);

% download: https://drive.google.com/drive/folders/1EWlkqu3qofB5kB5QLcPkrFyJxFm2zBec?usp=sharing
DATA_DIRECTORY = fullfile(MAIN_DIRECTORY, 'data');
DTM_DIRECTORY  = fullfile(DATA_DIRECTORY, 'ExportDMR_BielyKriz');
PC_DIRECTORY   = fullfile(DATA_DIRECTORY, 'ExportMB_BielyKriz');

cd(MAIN_DIRECTORY);

addpath(MAIN_DIRECTORY);

colorMap = load('las_colormap.txt')/255; % color map for classified point cloud

COLOR_GROUND     = colorMap(2,:);
COLOR_VEG_LOW    = colorMap(3,:);
COLOR_VEG_MEDIUM = colorMap(4,:);
COLOR_VEG_HIGH   = colorMap(5,:);

%% DATA READING - DMR
cd(DTM_DIRECTORY);

[DTM, rasterReference] = readgeoraster('dmr.tif', 'OutputType', 'double');
infoDMR = georasterinfo('dmr.tif');

DTM = standardizeMissing(DTM, infoDMR.MissingDataIndicator);
DTM = flipud(DTM);
rasterReference.ColumnsStartFrom = 'south';

%% PLOT - DMR
figure('Name', 'Digital Terrain Model data')
title 'DTM data'
mapshow(DTM, rasterReference, 'DisplayType', 'surface');
demcmap(DTM, 256)
view(3);
daspect([1 1 1]) % cim je posledne cislo blizsie ku 0, tym viac sa skaluje "z" suradnica
axis normal

%% DATA READING - PC
cd(PC_DIRECTORY);
lasFiles = dir('*.las'); % find all *.las files
lasCount = length(lasFiles); % number of found *.las files

ptCloud = cell(lasCount, 1);
ptAttributes = cell(lasCount, 1);
colorData = cell(lasCount, 1);

nPoints = 0;

for i = 1:lasCount
    lasReader = lasFileReader(lasFiles(i).name);
    [ptCloud{i}, ptAttributes{i}] = readPointCloud(lasReader, 'Attributes', 'Classification');
    colorData{i} = reshape(label2rgb(ptAttributes{i}.Classification, colorMap, 'k'), [], 3);
	nPoints = nPoints + ptCloud{i}.Count;
end

% average number of points
fprintf("Average point count per pixel: %.2f\n", nPoints / nnz(~isnan(DTM(:))));

%% PLOT - PC
% 2 -> Ground
% 3 -> Low Vegetation
% 4 -> Medium Vegetation
% 5 -> High Vegetation
selectedClasses = [3, 4, 5];
	
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

%% Data preprocessor
cd(MAIN_DIRECTORY);

preprocessor = dataPreprocessor(DTM, rasterReference, ptCloud, ptAttributes);
preprocessor = preprocessor.meshPlane();
preprocessor = preprocessor.filterPointCloud();
preprocessor = preprocessor.normalizePtCloud("method","DTM");

%% PLOT ORIGINAL AND NORMALIZED POINT CLOUD
figure("Name","Original Point Cloud")
preprocessor.plotPtCloud3D(colorMap, [2, 3, 4, 5],"useData","original")
%%
figure("Name","Normalized Point Cloud")
preprocessor.plotPtCloud3D(colorMap, [2, 3, 4, 5],"useData","normalized")
% handle.plotPtCloud3D(colorMap, [5],"useData","normalized")

%% COMPUTE ALL METRICS FOR GIVEN POINT CLOUD
featureExtractor = dataFeatureExtractor(rasterReference,... % original RasterReference
										preprocessor.ptCloud_norm,... % normalized point cloud
										preprocessor.ptAttributes_norm,... % and its attributes
										10); % desired pixel size in meters

%%
featureExtractor = featureExtractor.meshPlane();
[featureExtractor, time_1] = featureExtractor.computeMetricRasters();
% % % % [featureExtractor, time_2] = featureExtractor.computeMetricRastersParallel(6);

%% Kopia povodnych metrik
metricRastersCopy = featureExtractor.metricsRasters;

%% PLOT SELECTED METRIC
figure
featureExtractor.plotMetricRaster("plotData","Hmedian")
colormap jet
% colormap hot
colorbar
axis equal
axis xy

%% PLOT ALL METRICS
featureExtractor.plotAllMetricRasters("colormap","jet");

%% CLIP RASTER VALUES IF NECESSARY
% clip rasters: Coeff_var_z, Hkurt, Hskew, niekdy mozno aj Hvar
featureExtractor = featureExtractor.clipMetricRaster("clipData","Coeff_var_z","percentile",97.5);
featureExtractor = featureExtractor.clipMetricRaster("clipData","Hkurt","percentile",97.5);
featureExtractor = featureExtractor.clipMetricRaster("clipData","Hskew","percentile",97.5);
featureExtractor = featureExtractor.clipMetricRaster("clipData","Hvar","percentile",97.5);

%% skusanie vykreslenia
figure
imagesc([featureExtractor.x1 featureExtractor.x2],...
		[featureExtractor.y1 featureExtractor.y2], temp2,...
	'AlphaData', featureExtractor.alphaData)
title '99th perc'
colormap jet
colorbar
axis equal
axis xy

%% EXPORT SELECTED METRIC TO .TIF FILE
featureExtractor.exportMetricRaster("data1","exportLayer","Coeff_var_z");

%% EXPORT ALL METRICS TO .TIF FILES
featureExtractor.exportAllMetricRasters("data_BielyKriz_10x10m");
%%
% [H, R] = readgeoraster('Hmax.tif', 'OutputType', 'double');

% H = flipud(H);

%% Export color palette
% c = jet;
% c = c*255;
% writematrix(c, "jet_palette.txt","Delimiter",' ');
% 



