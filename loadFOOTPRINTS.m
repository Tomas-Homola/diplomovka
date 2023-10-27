function FOOTPRINTS_curves = loadFOOTPRINTS(SHP_filePath)
	% Function that read .SHP file with curves for all LOTS, then stores
	% the curves as polyshape objects and returns them in cell LOT_curves
	% with LOT number
	arguments
		SHP_filePath string
	end

	shapeFile = shaperead(SHP_filePath);

	numCurves = numel(shapeFile);
	FOOTPRINTS_curves = cell(numCurves, 2);

	for i = 1:numCurves
		% get polygon coordinates
		X = shapeFile(i).X;
		Y = shapeFile(i).Y;

		% get file name
		fileName = shapeFile(i).filename;
		fileName = strcat(fileName(1:end-4),'_jtsk03_bpv.laz');
		
		% save polygon and LOT number to cell
		FOOTPRINTS_curves{i,1} = polyshape(X, Y);
		FOOTPRINTS_curves{i,2} = fileName;
	end
end



