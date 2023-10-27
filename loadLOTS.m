function LOT_curves = loadLOTS(SHP_filePath)
	% Function that read .SHP file with curves for all LOTS, then stores
	% the curves as polyshape objects and returns them in cell LOT_curves
	% with LOT number
	arguments
		SHP_filePath string
	end

	shapeFile = shaperead(SHP_filePath);

	numCurves = numel(shapeFile);
	LOT_curves = cell(numCurves, 2);

	for i = 1:numCurves
		% get polygon coordinates
		X = shapeFile(i).X;
		Y = shapeFile(i).Y;

		% get LOT number
		lotNum = str2num(shapeFile(i).Cislo_loka);
		
		% save polygon and LOT number to cell
		LOT_curves{lotNum,1} = polyshape(X, Y);
		LOT_curves{lotNum,2} = shapeFile(i).Cislo_loka;
	end
end

