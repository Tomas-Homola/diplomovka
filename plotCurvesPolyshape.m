function plotCurvesPolyshape(curves)
	arguments
		curves cell
	end

	if (size(curves,2) == 2)
		hold on
		for i = 1:length(curves)
			plot(curves{i,1});
		end
		hold off
	else
		hold on
		for i = 1:length(curves)
			plot(curves{i});
		end
		hold off
	end

end

