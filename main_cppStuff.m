close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);

% download: https://drive.google.com/drive/folders/1EWlkqu3qofB5kB5QLcPkrFyJxFm2zBec?usp=sharing
DATA_DIRECTORY = fullfile(MAIN_DIRECTORY, 'data');
CURVES_DIRECTORY = fullfile(DATA_DIRECTORY, 'curves');

SHAPEFILES_DIRECTORY = fullfile(DATA_DIRECTORY,'_shapeFiles');
SHP_ALL_LOTS = fullfile(SHAPEFILES_DIRECTORY,...
						'prehlad_lokalit_lls_1_cyklus/prehlad_lokalit_lls_1_cyklus.shp');

SHP_FOOTPRINTS = cell(42, 1);
MAT_FOOTPRINTS = cell(42, 1);
lot_i = "";
for i = 1:42
	if (i <= 9)
		lot_i = "LOT0" + num2str(i);
	else
		lot_i = "LOT" + num2str(i);
	end
	
	SHP_FOOTPRINTS{i} = fullfile(SHAPEFILES_DIRECTORY,...
							  "footprints_s-jtsk03_shp/" + lot_i + "/" + lot_i + "_las.shp");

	MAT_FOOTPRINTS{i} = fullfile(SHAPEFILES_DIRECTORY,"footprints_mat/" + num2str(i) + ".mat");

end

LOT_DIR = fullfile(MAIN_DIRECTORY, "data/lazFiles");

cd(MAIN_DIRECTORY);

addpath(MAIN_DIRECTORY);

lasColorMap = load('las_colormap.txt')/255; % color map for classified point cloud

COLOR_GROUND     = lasColorMap(2,:);
COLOR_VEG_LOW    = lasColorMap(3,:);
COLOR_VEG_MEDIUM = lasColorMap(4,:);
COLOR_VEG_HIGH   = lasColorMap(5,:);

%%
% LOAD .KML CURVES
[KMLFiles, KMLFilePath] = uigetfile('*.kml', 'Select KML file', CURVES_DIRECTORY,'MultiSelect','off');
KMLFiles = string(KMLFiles);
KMLFilePath = string(KMLFilePath);
tic
curves = cell(numel(KMLFiles), 1);

for i = 1:numel(curves)
	coords = readKML(KMLFilePath + KMLFiles(i));
	if isempty(coords)
		continue;
	end
	lon = coords(:,1)';
	lat = coords(:,2)';
	h = coords(:,3)';
	
	[x, y] = gps_to_JTSK03_transformation(lat, lon, h);
	curves{i}.poly = polyshape(x, y);
	curves{i}.name = KMLFiles(i);

	fprintf("Curve %d/%d loaded\n", i, numel(curves));
end
lazTime = toc;
clear h x y lat lon coords

fprintf("All KML curves loaded in %.2f s.\n", lazTime);

cd(MAIN_DIRECTORY);
h = 10; n = 1;

dataPP = dataPreprocessorCurve(curves{1}.poly, h, n);

% FIND NAMES OF .LAZ FILES FOR EACH LOADED CURVE
foundLots = []; % cell of strings of found LOTS where the curves{i} lies inside
foundLazFiles = cell(length(curves), 1);
LOTS_curves = loadLOTS(SHP_ALL_LOTS); % table of curves for
fprintf("LOTs curves loaded\n");

tic
for i = 1:length(curves) % iterate through all .KML curves
	foundLots = []; % cell of strings of found LOTS where the curves{i} lies inside
	% find indices of LOTs where the i-th .KML curve belongs to
	for j = 1:length(LOTS_curves) % iterate through all LOT curves
		% find intersection of i-th .KML curve with j-th LOT
		intersection = intersect(dataPP.Omega, LOTS_curves{j}.poly);

		if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
			continue;
		end

		foundLots = [foundLots, str2num(LOTS_curves{j}.lotNum)];

	end % end j
	
	% find names of .LAZ files to load for the i-th curve
	foundLazFiles{i} = {};
	lasFileFound = 0;
	FOOTPRINTS_curves = [];
	
	if isempty(foundLots)
		fprintf("Curve %s outside SVK\n", curves{i}.name);
		continue;
	end
	
	for j = 1:length(foundLots) % iterate over all found LOTs
% 		FOOTPRINTS_curves = loadFOOTPRINTS(SHP_FOOTPRINTS{foundLots(j)});
		load(MAT_FOOTPRINTS{foundLots(j)});
% 		fprintf("Footprints loaded\n");
		
		for k = 1:length(FOOTPRINTS_curves) % iterate over all footprints within found j-th LOT

			intersection = intersect(dataPP.Omega, FOOTPRINTS_curves{k}.poly);

			if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
				continue;
			end
		
			% save found .LAZ file name
			lasFileFound = lasFileFound + 1;
			foundLazFiles{i}{lasFileFound}.poly = FOOTPRINTS_curves{k}.poly;
			foundLazFiles{i}{lasFileFound}.lazName = FOOTPRINTS_curves{k}.fileName;
			foundLazFiles{i}{lasFileFound}.lazPath = strrep(fullfile(LOT_DIR, FOOTPRINTS_curves{k}.fileName),'\', '/');
		end
		
	end

	fprintf("Curve %d/%d done\n", i, length(curves));
% 	fprintf("Found .LAZ files:\n");
% 	for j = 1:length(foundLazFiles{i})
% 		fprintf("%s\n", foundLazFiles{i}{j});
% 	end
end % end i

lazTime = toc;
fprintf("LAZ files for all curves found, elapsed time: %.2f s\n", lazTime);

% PRINT FOUND .LAZ FILES
for i = 1:length(foundLazFiles)
	fprintf("\nFound .LAZ files for curve '%s':\n", curves{i}.name);
	for j = 1:length(foundLazFiles{i})
	fprintf("%s\n", foundLazFiles{i}{j}.lazName);
	end
end

%% PLOT KML CURVE AND FOUND FOOTPRINTS
% plot(curves{2},"LineWidth",5,"EdgeColor","red")
c = 1;
figure
axis equal
hold on
plot(dataPP.Omega,"LineWidth",4,"EdgeColor","magenta")
for i = 1:length(foundLazFiles{c})
	plot(foundLazFiles{c}{i}.poly,"FaceAlpha",0.05);
	pause(0.5)
end
plot(curves{c}.poly,"LineWidth",4,"EdgeColor","red","FaceAlpha",0)
hold off

%% Create out file for C++
outFile = fopen(curves{1}.name + "_data.txt", 'w');

fprintf(outFile, "%d;%d\n", dataPP.nx, dataPP.ny);
fprintf(outFile, "%d;%d\n", dataPP.nx/h, dataPP.ny/h);
fprintf(outFile, "%.10f;%.10f\n", dataPP.x1, dataPP.y2);
fprintf(outFile, "%d\n", length(foundLazFiles{c}));

for i = 1:length(foundLazFiles{c}) % nacitanie PC cez .MAT subory pre krivku "c"
	fprintf(outFile, "%s\n", foundLazFiles{c}{i}.lazPath);
end

fclose(outFile);

%% Visualize computed metrics
cd(DATA_DIRECTORY);

filename = '2024_02_25_18-14-35_91E0_bodikyprihradzi1_final.kml_data_h=10m.tif';
[metrics, rasterReference] = readgeoraster(filename, 'OutputType', 'double');
infoMetrics = georasterinfo(filename);

metrics = standardizeMissing(metrics, infoMetrics.MissingDataIndicator);
metrics = flipud(metrics);
rasterReference.ColumnsStartFrom = 'south';

[X, Y] = rasterReference.worldGrid();

band = 5;
figure('Name', 'Metrics Bodiky')
title 'Bodiky data'
alphaData = ones(size(X));
alphaData(isnan(metrics(:,:,band))) = 0;
imagesc(X(1,:), Y(:,1), metrics(:,:,band), 'AlphaData', alphaData)
set(gca,'YDir','normal') % spravne hodnoty na y osi
colormap bone
colorbar
axis equal
% plotCurvesPolyshape(curves)

%% Find all necessarry laz files for all curves
uniqueLaz = cell(1);
uniqueLaz{1} = '';
unique = 1;
for i = 1:length(foundLazFiles)
	for j = 1:length(foundLazFiles{i})
		if (~ismember(foundLazFiles{i}{j}.lazName, uniqueLaz))
			uniqueLaz{unique} = foundLazFiles{i}{j}.lazName;
			unique = unique + 1;
		end

	end
end

uniqueLaz = uniqueLaz';

%% Find missing laz files inside laz files folder
cd(LOT_DIR)

missingFiles = cell(1);
missingFiles{1} = '';
found = 1;
for i = 1:length(uniqueLaz)
	if (~isfile(uniqueLaz{i})) % if laz not found
		missingFiles{found} = uniqueLaz{i};
		found = found + 1;
	end
end

missingFiles = missingFiles'

%% Copy missing files - current folder path to hard-drive
for i = 1:length(missingFiles)
	fprintf("file %d/%d\n", i, length(missingFiles));
	if (isfile(missingFiles{i}))
		if (~isfile(fullfile(MAIN_DIRECTORY, "data/lazFiles", missingFiles{i})))
			copyfile(missingFiles{i}, fullfile(MAIN_DIRECTORY, "data/lazFiles", missingFiles{i}))
		end
	else
		fprintf("file %d not found\n", i);
	end
end

fprintf("Copy files done\n");


