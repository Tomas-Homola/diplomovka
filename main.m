% close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);

% download: https://drive.google.com/drive/folders/1EWlkqu3qofB5kB5QLcPkrFyJxFm2zBec?usp=sharing
DATA_DIRECTORY = fullfile(MAIN_DIRECTORY, 'data');
DTM_DIRECTORY  = fullfile(DATA_DIRECTORY, 'ExportDMR_porovnanie');
PC_DIRECTORY   = fullfile(DATA_DIRECTORY, 'ExportMB_porovnanie');
CURVES_DIRECTORY = fullfile(DATA_DIRECTORY, 'curves');

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
selectedClasses = [2 5];
	
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
preprocessor.plotPtCloud3D(colorMap, [2, 3, 5],"useData","original")
%%
figure("Name","Normalized Point Cloud")
preprocessor.plotPtCloud3D(colorMap, [2, 3, 5],"useData","normalized")
% handle.plotPtCloud3D(colorMap, [5],"useData","normalized")


%% COMPUTE ALL METRICS FOR GIVEN POINT CLOUD
featureExtractor = dataFeatureExtractor(rasterReference,... % original RasterReference
										preprocessor.ptCloud_norm,... % normalized point cloud
										preprocessor.ptAttributes_norm,... % and its attributes
										5); % desired pixel size in meters

featureExtractor = featureExtractor.meshPlane();
%% 
[featureExtractor, time_1] = featureExtractor.computeMetricRasters();
% % % % [featureExtractor, time_2] = featureExtractor.computeMetricRastersParallel(6);

%% Vykreslenie rezu PC
% H = flipud(DTM);
H = flipud(featureExtractor.metricsRasters.Hskew);
figure
Alpha = ones(size(H));
Alpha(isnan(H)) = 0;
imagesc([preprocessor.xc(1) preprocessor.xc(end)], [preprocessor.yc(end) preprocessor.yc(1)], ...
		H, 'AlphaData', Alpha)
set(gca,'YDir','normal') % spravne hodnoty na y osi
% demcmap(H, 256);
colormap jet
colorbar
A = ginput(1);
hold on
plot(A(1), A(2), '*k','MarkerSize',15)
B = ginput(1);
plot([A(1), B(1)], [A(2), B(2)], '*-k', 'LineWidth', 4, 'MarkerSize',15)
drawnow
pause(0.5)
hold off

k = (B(2) - A(2))/(B(1) - A(1));
q = A(2) - k * A(1);

% fprintf("k: %.5f\nq: %.5f\n", k, q)
d = @(X) abs(-k*X(:,1) + X(:,2) - q)/sqrt(k*k + 1);

width = 5;
selectedClasses = [2 3 4 5 6 9];

figure('Name', 'Point Cloud original')
title 'Point Cloud original'
hold on
for i = 1:lasCount
	% original
	xy = ptCloud{i}.Location(:,1:2);
	t = (xy - A) / (B - A);
	selected = d(xy) < width & t <= 1 & t >= 0;
	
	classMember = ismember(ptAttributes{i}.Classification, selectedClasses);
		if any(classMember & selected)
			pcshow(ptCloud{i}.Location(classMember & selected, :), colorData{i}(classMember & selected, :), ...
				'MarkerSize', 20)
		end
end
hold off

PC_n = preprocessor.ptCloud_norm;
att_n = preprocessor.ptAttributes_norm';

figure('Name', 'Point Cloud normalized')
title 'Point Cloud normalized'
hold on
for i = 1:lasCount
	% normalized
	xy = PC_n{i}.Location(:,1:2);
	t = (xy - A) / (B - A);
	selected = d(xy) < width & t <= 1 & t >= 0;
	
	classMember = ismember(att_n{i}.Classification, selectedClasses);
		if any(classMember & selected)
			colorData_i = reshape(label2rgb(att_n{i}.Classification, colorMap, 'k'), [], 3);
			pcshow(PC_n{i}.Location(classMember & selected, :), colorData_i(classMember & selected, :), ...
				'MarkerSize', 20)
		end
end
hold off

%% Kopia povodnych metrik
metricRastersCopy = featureExtractor.metricsRasters;
alphaDataCopy = featureExtractor.alphaData;


%% LOAD CURVES
[KMLFiles, KMLFilePath] = uigetfile('*.kml', 'Select KML file', CURVES_DIRECTORY,'MultiSelect','on');
KMLFiles = string(KMLFiles);
KMLFilePath = string(KMLFilePath);

curves = cell(numel(KMLFiles), 1);

for i = 1:numel(curves)
	coords = readKML(KMLFilePath + KMLFiles(i));
	lon = coords(:,1)';
	lat = coords(:,2)';
	h = coords(:,3)';
	
	[x, y] = gps_to_JTSK03_transformation(lat, lon, h);
	curves{i} = curve();
	curves{i} = curves{i}.set_coordinates([x', y']);
end

clear h x y lat lon coords

fprintf("Curves loaded.\n");

%% PLOT SELECTED METRIC
figure('WindowState','maximized')
featureExtractor.plotMetricRaster("plotData","Hmax")
plotCurves(curves,"Color",'magenta','MarkerSize',20)
colormap jet
% colormap hot
colorbar
axis equal
axis xy

%% PLOT ALL METRICS
% pre ziadne krivky je curves = {}
featureExtractor.plotAllMetricRasters("colormap","jet","plotCurves",{});

%% CLIP RASTER VALUES IF NECESSARY
% clip rasters: Coeff_var_z, Hkurt, Hskew, niekdy mozno aj Hvar
featureExtractor = featureExtractor.clipMetricRaster("clipData","Coeff_var_z","percentile",97.5);
featureExtractor = featureExtractor.clipMetricRaster("clipData","Hkurt","percentile",97.5);
featureExtractor = featureExtractor.clipMetricRaster("clipData","Hskew","percentile",97.5);
featureExtractor = featureExtractor.clipMetricRaster("clipData","Hvar","percentile",97.5);
featureExtractor = featureExtractor.clipMetricRaster("clipData","DAM_z","percentile",97.5);

%% EXPORT SELECTED METRIC TO .TIF FILE
featureExtractor.exportMetricRaster("data1","exportLayer","Coeff_var_z");

%% EXPORT ALL METRICS TO .TIF FILES
fileName = "data_BielyKriz_5x5m";
featureExtractor.exportAllMetricRasters(fileName);

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

%%
% figure
% imagesc([preprocessor.x1 preprocessor.x2],...
% 		[preprocessor.y1 preprocessor.y2], preprocessor.DTM)
% title '99th perc'
% plotCurves(curves,"Color",'magenta')
% colormap jet
% colorbar
% axis equal
% axis xy


%%
% [H, R] = readgeoraster('Hmax.tif', 'OutputType', 'double');

% H = flipud(H);

%% Export color palette
% c = jet;
% c = c*255;
% writematrix(c, "jet_palette.txt","Delimiter",' ');
% 

%% NACITANIE .SHP SUBORU
% [shpFileName, shpFilePath] = uigetfile('*.shp', 'Select SHP file', DATA_DIRECTORY); % .osm file path
% shpFilePath = fullfile(shpFilePath, shpFileName);
%%
% shapeFileInfo = shapeinfo(shpFilePath);
% shapeFile = shaperead(shpFilePath);

%%
% BBcoords = shapeFileInfo.BoundingBox;
% X = [BBcoords(1,1) BBcoords(2,1) BBcoords(2,1) BBcoords(1,1)];
% Y = [BBcoords(1,2) BBcoords(1,2) BBcoords(2,2) BBcoords(2,2)];
% boundingBox = polyshape(X, Y);
% 
% numCurves = numel(shapeFile);
% curves = cell(numCurves, 1);
% 
% figure
% plot(boundingBox,"FaceAlpha",0);
% hold on
% for i = 1:numCurves
% 	X = shapeFile(i).X;
% 	Y = shapeFile(i).Y;
% 	curves{i} = polyshape(X, Y);
% 	
% 	plot(curves{i});
% 	pause(0.05);
% end
% 
% hold off


