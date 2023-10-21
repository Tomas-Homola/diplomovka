function plotCurves(curvesCell, options)
	arguments
		curvesCell (:,1) cell
		options.Color = 'red'
		options.MarkerSize = 10;
	end

	hold on
	for i = 1:numel(curvesCell)
		curvesCell{i}.plot("Color", options.Color, "MarkerSize", options.MarkerSize)
	end
	hold off
end

