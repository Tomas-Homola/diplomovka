function coords = readKML(fileName)
	if exist(fileName, 'file') ~= 2
		error("File does not exist");
	end

	% open file
	fileID = fopen(fileName, 'r');

	% read all lines from the kml file
	kmlLines = textscan(fileID, '%s', 'Delimiter','\n');
	kmlLines = kmlLines{1};

	% close file
	fclose(fileID);

	% find where the coordinates start and end
	coordsStart = find(contains(kmlLines, '<coordinates>'), 1, "first");
	coordsEnd   = find(contains(kmlLines, '</coordinates>'), 1, "first");
	
	% save lines containing coordinates
	coordsLines = kmlLines(coordsStart+1:coordsEnd-1);
	
	% allocate space for coordinates
	coords = zeros(numel(coordsLines), 3);

	% read coordinates and save them
	for i = 1:numel(coordsLines)
		values = str2num(coordsLines{i});
		
		coords(i,:) = values;
	end

end









