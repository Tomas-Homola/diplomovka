function plotCurvesPolyshape(curves, options)
arguments
	curves cell
	options.FaceAlpha = 1
	options.Incr = 1
end

hold on
for i = 1:options.Incr:length(curves)
	plot(curves{i}.poly,"FaceAlpha", options.FaceAlpha);
end
hold off

end

