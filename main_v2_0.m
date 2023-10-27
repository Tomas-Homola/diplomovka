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
lot_i = "";
for i = 1:42
	if (i <= 9)
		lot_i = "LOT0" + num2str(i);
	else
		lot_i = "LOT" + num2str(i);
	end
	
	SHP_FOOTPRINTS{i} = fullfile(SHAPEFILES_DIRECTORY,...
							  "footprints_s-jtsk03_shp/" + lot_i + "/" + lot_i + "_las.shp");

end

cd(MAIN_DIRECTORY);

addpath(MAIN_DIRECTORY);

colorMap = load('las_colormap.txt')/255; % color map for classified point cloud

COLOR_GROUND     = colorMap(2,:);
COLOR_VEG_LOW    = colorMap(3,:);
COLOR_VEG_MEDIUM = colorMap(4,:);
COLOR_VEG_HIGH   = colorMap(5,:);

%% LOAD .KML CURVES
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
	curves{i} = polyshape(x, y);
end

clear h x y lat lon coords

fprintf(".KML curves loaded.\n");

% plot loaded curves
figure
plotCurvesPolyshape(curves)
axis equal

%% FIND NAMES OF .LAZ FILES FOR EACH LOADED CURVE
foundLots = []; % cell of strings of found LOTS where the curves{i} lies inside
foundLazFiles = cell(length(curves), 1);
LOTS_curves = loadLOTS(SHP_ALL_LOTS); % table of curves for
fprintf("LOTs curves loaded\n");

for i = 1:length(curves) % iterate through all .KML curves

	% find indices of LOTs where the i-th .KML curve belongs to
	for j = 1:length(LOTS_curves) % iterate through all LOT curves
		% find intersection of i-th .KML curve with j-th LOT
		intersection = intersect(curves{i}, LOTS_curves{j,1});

		if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
			continue;
		end

		foundLots = [foundLots; str2num(LOTS_curves{j,2})];

	end % end j
	
	fprintf("aaaa\n");
	% find names of .LAZ files to load for the i-th curve
	for j = 1:length(foundLots)
		FOOTPRINTS_curves = loadFOOTPRINTS(SHP_FOOTPRINTS{foundLots(j)});
		
		% TODO: zistenie, kam patri i-ta krivka vramci footprints
	end


end % end i








