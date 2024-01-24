close all;
clear;
clc;

%% DIRECTORY STUFF
[MAIN_DIRECTORY,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(MAIN_DIRECTORY);
cd("data\lazFiles\");

lazFiles = dir("*.laz");

lazNames = cell(0);
matNames = cell(0);
for i = 1:length(lazFiles)
	lazNames{i} = lazFiles(i).name;
	matNames{i} = strcat( lazNames{i}(1:end-3), 'mat' );
end
lazNames = lazNames';
matNames = matNames';

ptCloud_i = []; 
ptAttributes_i = [];

for i = 1:length(lazNames)
	% ak este neexistuje dany .mat subor, tak sa vytvori
	if (~isfile(matNames{i}))
		fprintf("Started file %d/%d... ", i, length(lazNames))
		tic;
		lasReader = lasFileReader(lazNames{i});
		[ptCloud_i, ptAttributes_i] = readPointCloud(lasReader,...
									'Attributes', 'Classification');

		save(matNames{i}, "ptCloud_i", "ptAttributes_i");
		pause(2);
		ptCloud_i = []; %#ok
		ptAttributes_i = []; %#ok

		time = toc;
		fprintf("done in %.3f\n", time);
	end
end

fprintf("Converting .laz to .mat files finished\n");

