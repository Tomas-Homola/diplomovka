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
tic
curves = cell(numel(KMLFiles), 1);

for i = 1:numel(curves)
	coords = readKML(KMLFilePath + KMLFiles(i));
	lon = coords(:,1)';
	lat = coords(:,2)';
	h = coords(:,3)';
	
	[x, y] = gps_to_JTSK03_transformation(lat, lon, h);
	curves{i}.poly = polyshape(x, y);
	curves{i}.name = KMLFiles(i);

	fprintf("Curve %d/%d loaded\n", i, numel(curves));
end
time = toc;
clear h x y lat lon coords

fprintf("All KML curves loaded in %.2f s.\n", time);
%%
% plot loaded curves
figure
plotCurvesPolyshape(curves)
axis equal

%% Pomocne vykreslenie footprints kriviek
% figure
% hold on
% for i = 1:length(FOOTPRINTS_curves)
% 	plot(FOOTPRINTS_curves{i}.poly)
% end
% hold off
% axis equal

%% FIND NAMES OF .LAZ FILES FOR EACH LOADED CURVE
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
		intersection = intersect(curves{i}.poly, LOTS_curves{j}.poly);

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

			intersection = intersect(curves{i}.poly, FOOTPRINTS_curves{k}.poly);

			if isempty(intersection.Vertices) % if i-th .KML curve does not intersect with j-th LOT curve -> continue
				continue;
			end
		
			% save found .LAZ file name
			lasFileFound = lasFileFound + 1;
			foundLazFiles{i}{lasFileFound}.poly = FOOTPRINTS_curves{k}.poly;
			foundLazFiles{i}{lasFileFound}.lazName = FOOTPRINTS_curves{k}.fileName;
		end
		
	end

	fprintf("Curve %d/%d done\n", i, length(curves));
% 	fprintf("Found .LAZ files:\n");
% 	for j = 1:length(foundLazFiles{i})
% 		fprintf("%s\n", foundLazFiles{i}{j});
% 	end
end % end i

time = toc;
fprintf("LAZ files for all curves found, elapsed time: %.2f s\n", time);

%% PRINT FOUND .LAZ FILES
for i = 1:4
	fprintf("\nFound .LAZ files for curve '%s':\n", curves{i}.name);
	for j = 1:length(foundLazFiles{i})
	fprintf("%s\n", foundLazFiles{i}{j}.lazName);
	end
end

%% PLOT KML CURVE AND FOUND FOOTPRINTS
% plot(curves{2},"LineWidth",5,"EdgeColor","red")
c = 10;
figure
axis equal
hold on
plot(curves{c}.poly,"LineWidth",5,"EdgeColor","red")
for i = 1:length(foundLazFiles{c})
	plot(foundLazFiles{c}{i}.poly,"FaceAlpha",0.05);
	pause(1)
end
hold off

%%
[xlim, ylim] = curves{1}.poly.boundingbox;
BB_width = abs(xlim(2)-xlim(1)); % sirka bounding boxu kml krivky
BB_height = abs(ylim(2) - ylim(1)); % vyska bounding boxu kml krivky

h = 10; % pozadovana velkost pixelu v metroch
n = 1; % kolko nasobok "h" pridat ku sirke/vyske \Omega

% vypocet sirky a vysky pre \Omega
omega_width = ceil(BB_width / h) * h + n * h;
omega_height = ceil(BB_height / h) * h + n * h;

% fprintf("BB width: %.2f -> omega width: %.2f\nBB height: %.2f -> omega height: %.2f\n",...
% 	BB_width, omega_width, BB_height, omega_height);

% offset v smeroch x a y -> sluzi na vypocet suradnic pre vrcholy \Omega
offset_x = (omega_width - BB_width) / 2;
offset_y = (omega_height - BB_height) / 2;

% vypocet suradnic pre vrcholy \Omega
x_new = [xlim(1) - offset_x, xlim(2) + offset_x, xlim(2) + offset_x, xlim(1) - offset_x];
y_new = [ylim(2) + offset_y, ylim(2) + offset_y, ylim(1) - offset_y, ylim(1) - offset_y];

omega = polyshape(x_new, y_new);
[xlim_new, ylim_new] = omega.boundingbox;

% ulozenie suradnic kml krivky -> kvoli funkcii inpolygon
x = curves{1}.poly.Vertices(:,1);
y = curves{1}.poly.Vertices(:,2);

% randX = xlim_new(1) + (xlim_new(2) - xlim_new(1)) * rand(5000,1);
% randY = ylim_new(1) + (ylim_new(2) - ylim_new(1)) * rand(5000,1);
% randPoints = [randX, randY];
% 
% [in, on] = inpolygon(randX, randY, x, y);

figure
hold on
% scatter(randX(in),randY(in),'.g')
% scatter(randX(~in),randY(~in),0.1,'.r')
% scatter(randX(on),randY(on),20,'*b')

plot(curves{1}.poly, "EdgeColor", "blue", "FaceAlpha", 0.05, "FaceColor","blue")
plot(omega,"FaceAlpha",0,"LineStyle","-","LineWidth",2)
axis equal

% vytvorenie vypoctovej oblasti
x1 = xlim_new(1);
x2 = xlim_new(2);
y1 = ylim_new(1);
y2 = ylim_new(2);

nx = omega_width / h;
ny = omega_height / h;

xc = linspace(x1 + h/2, x2 - h/2, nx);
yc = linspace(y1 + h/2, y2 - h/2, ny);
[Xc, Yc] = meshgrid(xc, yc);

% pomocne premenne pre vykreslenie mriezky
X = x1:h:x2;
Y = y1:h:y2;

for i = 1:length(X)
	plot([X(i) X(i)],[Y(1) Y(end)],...
		'Color', "#EDB120", 'LineWidth', 0.5,'LineStyle','-') %y grid lines
end

for i = 1:length(Y)
	plot([X(1) X(end)],[Y(i) Y(i)],...
		'Color', "#EDB120", 'LineWidth', 0.5,'LineStyle','-') %x grid lines
end

% najdenie tych pixelov, ktorych stredy lezia v kml krivke
[XYin, XYon] = inpolygon(Xc, Yc, x, y);

scatter(Xc(XYin), Yc(XYin),5, "*g")
scatter(Xc(~XYin), Yc(~XYin),5, "*r")

hold off

legend('KML curve','\Omega', '\Omega discr')

fprintf("In square for n = %.2f -> %d pixels\n", n, nnz(XYin))




