classdef dataFeatureExtractorCurve
	% Class for feature extraction from point cloud data using some
	% reference curve
	% 1) load normalized point cloud from dataPreprocessor and RasterReference
	% 2) create MESH with pixel size h x h meters
	% 3) find indices within the new mesh
	% 4) compute metrics and save into preallocated space
	
	properties
		RR map.rasterref.MapCellsReference % RasterReference of the new mesh

		curve polyshape					% curve saved as polyshape object
		Omega polyshape

		% 2D mesh parameters
		x1, x2, y1, y2                  % x/y pixel centers range
        nx, ny                          % number of pixels in x/y direction
        h                               % desired pixel size
		n
        x, y                            % pixel edge coordinates
        xc, yc, Xc, Yc                  % pixel center coordinates (and meshgrids)
		nPixels							% number of mesh pixels
		isPixelInCurve					% logical vector of mesh pixels within given curve
		alphaData

		% Point cloud stuff
		lasCount						% number of point clouds
		ptCloud						    % cell of normalized point cloud
		ptAttributes					% cell of classification for all points clouds
		ptCloud_pixels					% pixel indices of all points within the mesh

		% Important point cloud classes
		GROUND = 2
		VEG_LOW = 3
		VEG_MED = 4
		VEG_HIGH = 5

		% structure for rasters of individual metrics computed from point cloud
		metricsRasters
		representativeMetrics			% array for representative metrics from within given curve
		% 
		plotTitles
		exportTitles
	end
	
	methods
		function this = dataFeatureExtractorCurve(curve, h, n, ptCloud_norm, ptAttributes_norm)
			arguments
				curve polyshape
				h (1,1) {mustBeNumeric, mustBeInteger, mustBeGreaterThanOrEqual(h, 1)}
				n
				ptCloud_norm
				ptAttributes_norm
			end

			this.ptCloud = ptCloud_norm;
			this.ptAttributes = ptAttributes_norm;
			this.lasCount = length(ptCloud_norm);

			this.curve = curve;
			this.h = h;
			this.n = n;

			% MESH
			[xlim_c, ylim_c] = curve.boundingbox;
			BB_width  = abs(xlim_c(2) - xlim_c(1)); % sirka bounding boxu kml krivky
			BB_height = abs(ylim_c(2) - ylim_c(1)); % vyska bounding boxu kml krivky

			% vypocet sirky a vysky pre Omega
			Omega_width  = ceil(BB_width / h) * h + n * h;
			Omega_height = ceil(BB_height / h) * h + n * h;

			% offset v smeroch x a y -> sluzi na vypocet suradnic pre vrcholy Omega
			offset_x = (Omega_width - BB_width) / 2;
			offset_y = (Omega_height - BB_height) / 2;

			% vypocet suradnic pre vrcholy Omega
			x_new = [xlim_c(1) - offset_x, xlim_c(2) + offset_x, xlim_c(2) + offset_x, xlim_c(1) - offset_x];
			y_new = [ylim_c(2) + offset_y, ylim_c(2) + offset_y, ylim_c(1) - offset_y, ylim_c(1) - offset_y];

			this.Omega = polyshape(x_new, y_new);
			[xlim_o, ylim_o] = this.Omega.boundingbox;

			% vytvorenie vypoctovej oblasti
			this.x1 = xlim_o(1);
			this.x2 = xlim_o(2);
			this.y1 = ylim_o(1);
			this.y2 = ylim_o(2);

			this.nx = Omega_width / h;
			this.ny = Omega_height / h;
			this.nPixels = this.ny * this.nx;

			% create mesh
			this.xc = linspace(this.x1 + this.h/2, this.x2 - this.h/2, this.nx);
			this.yc = linspace(this.y1 + this.h/2, this.y2 - this.h/2, this.ny);
			[this.Xc, this.Yc] = meshgrid(this.xc, this.yc);

			% find pixels centers inside curve
			x_c = curve.Vertices(:,1);
			y_c = curve.Vertices(:,2);
			[this.isPixelInCurve, ~] = inpolygon(this.Xc(:), this.Yc(:), x_c, y_c);

			% auxilliaries for ploting mesh
			this.x = this.x1:this.h:this.x2;
			this.y = this.y1:this.h:this.y2;

			% create new RasterReference object for the new mesh
			this.RR = maprefcells([this.x1, this.x2], [this.y1, this.y2], [this.ny, this.nx],...
				"ColumnsStartFrom","north");
			this.RR.ProjectedCRS = projcrs(8353); % S-JTSK [JTSK03] Krovak East North
			% this.RR.worldGrid() -> tato funkcia dokaze vytvorit mesh aky treba
			
			this = setMetricRastersNames(this);

		end % end of constructor

		function this = setMetricRastersNames(this)
			% initialize structures for metric rasters and plot titles
			this.metricsRasters = struct;

			this.plotTitles = struct;
			this.plotTitles.Hmax        = 'Max of normalized height';
			this.plotTitles.Hmean       = 'Mean of normalized height';
			this.plotTitles.Hmedian     = 'Median of normalized height';
			this.plotTitles.Hp25        = '25th percentile of normalized height';
			this.plotTitles.Hp75	    = '75th percentile of normalized height';
			this.plotTitles.Hp95        = '95th percentile of normalized height';
			
			this.plotTitles.PPR         = 'Pulse penetration ratio';
			this.plotTitles.DAM_z		= 'Number of returns above mean height';
			this.plotTitles.BR_bellow_1 = 'Proportion of vegetation points bellow 1 m';
			this.plotTitles.BR_1_2      = 'Proportion of vegetation points between 1 m - 2 m';
			this.plotTitles.BR_2_3      = 'Proportion of vegetation points between 2 m - 3 m';
			this.plotTitles.BR_above_3  = 'Proportion of vegetation points above 3 m';
			this.plotTitles.BR_3_4      = 'Proportion of vegetation points between 3 m - 4 m';
			this.plotTitles.BR_4_5      = 'Proportion of vegetation points between 4 m - 5 m';
			this.plotTitles.BR_bellow_5 = 'Proportion of vegetation points bellow 5 m';
			this.plotTitles.BR_5_20     = 'Proportion of vegetation points between 5 m - 20 m';
			this.plotTitles.BR_above_20 = 'Proportion of vegetation points above 20 m';
			
			this.plotTitles.Coeff_var_z  = 'Coefficient of variation';
			this.plotTitles.Hkurt        = 'Kurtosis of normalized height';
			this.plotTitles.Hskew        = 'Skewness of normalized height';
			this.plotTitles.Hstd         = 'Standard deviation of normalized height';
			this.plotTitles.Hvar         = 'Variance of normalized height';

			% initialize structure for export titles
			this.exportTitles = struct;

			this.exportTitles.Hmax        = '1_Hmax_';
			this.exportTitles.Hmean       = '2_Hmean_';
			this.exportTitles.Hmedian     = '3_Hmedian_';
			this.exportTitles.Hp25        = '4_Hp25_';
			this.exportTitles.Hp75	      = '6_HP75_';
			this.exportTitles.Hp95        = '7_HP95_';
			
			this.exportTitles.PPR         = '8_PPR_';
			this.exportTitles.DAM_z		  = '9_Density_above_mean_z_';
			this.exportTitles.BR_bellow_1 = '10_BR_bellow_1_';
			this.exportTitles.BR_1_2      = '11_BR_1_2_';
			this.exportTitles.BR_2_3      = '12_BR_2_3_';
			this.exportTitles.BR_above_3  = '13_BR_above_3_';
			this.exportTitles.BR_3_4      = '14_BR_3_4_';
			this.exportTitles.BR_4_5      = '15_BR_4_5_';
			this.exportTitles.BR_bellow_5 = '16_BR_bellow_5_';
			this.exportTitles.BR_5_20     = '17_BR_5_20_';
			this.exportTitles.BR_above_20 = '18_BR_above_20_';
			
			this.exportTitles.Coeff_var_z  = '19_Coeff_var_z_';
			this.exportTitles.Hkurt       = '21_Hkurt_';
			this.exportTitles.Hskew       = '23_Hskew_';
			this.exportTitles.Hstd        = '24_Hstd_';
			this.exportTitles.Hvar        = '25_Hvar_';

		end

		function plotMesh(this, options)
			arguments
				this 
				options.Color = [0.3, 0.3, 0.3]
				options.LineWidth (1,1) {mustBeNumeric} = 4
			end
			
			hold on
			for i = 1:length(this.x) % oranzova = #EDB120
				plot([this.x(i) this.x(i)],[this.y(1) this.y(end)],...
					'Color', "black", 'LineWidth', 1,'LineStyle','-') %y grid lines
			end

			for i = 1:length(this.y)
				plot([this.x(1) this.x(end)],[this.y(i) this.y(i)],...
					'Color', "black", 'LineWidth', 1,'LineStyle','-') %x grid lines
			end

			plot(this.curve, "EdgeColor", "magenta", "FaceAlpha", 0, "LineWidth",2)
			plot(this.Omega,"FaceAlpha",0,"LineStyle","-","LineWidth",2)

			x_ = this.curve.Vertices(:,1);
			y_ = this.curve.Vertices(:,2);
			[XYin, ~] = inpolygon(this.Xc, this.Yc, x_ , y_);

			scatter(this.Xc(XYin),  this.Yc(XYin), 13, "og","filled", "MarkerEdgeColor","black")
			scatter(this.Xc(~XYin), this.Yc(~XYin),13, "or","filled", "MarkerEdgeColor","black")
			
			hold off
			
			
		end % end of function plotMesh

		function [this, timePassed] = computeMetricRasters(this, computeAllPixels)
			arguments
				this dataFeatureExtractorCurve
				computeAllPixels logical = 0;
			end
			% allocate space for metric rasters
			% ECOSYSTEM HEIGHT
			this.metricsRasters.Hmax = NaN(this.ny, this.nx);
			this.metricsRasters.Hmean = NaN(this.ny, this.nx);
			this.metricsRasters.Hmedian = NaN(this.ny, this.nx);
			this.metricsRasters.Hp25 = NaN(this.ny, this.nx);
			this.metricsRasters.Hp75 = NaN(this.ny, this.nx);
			this.metricsRasters.Hp95 = NaN(this.ny, this.nx);

			% ECOSYSTEM COVER
			this.metricsRasters.PPR = NaN(this.ny, this.nx);
			this.metricsRasters.DAM_z = NaN(this.ny, this.nx);
			this.metricsRasters.BR_bellow_1 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_1_2 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_2_3 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_above_3 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_3_4 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_4_5 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_bellow_5 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_5_20 = NaN(this.ny, this.nx);
			this.metricsRasters.BR_above_20 = NaN(this.ny, this.nx);

			% ECOSYSTEM STRUCTURAL COMPLEXITY
			this.metricsRasters.Coeff_var_z = NaN(this.ny, this.nx);
			this.metricsRasters.Hkurt = NaN(this.ny, this.nx);
			this.metricsRasters.Hskew = NaN(this.ny, this.nx);
			this.metricsRasters.Hstd = NaN(this.ny, this.nx);
			this.metricsRasters.Hvar = NaN(this.ny, this.nx);
				
			time = tic;

			% find pixel indices for all points within mesh
			this.ptCloud_pixels = cell(this.lasCount, 1);
			for i = 1:this.lasCount
				xpt = this.ptCloud{i}.Location(:, 1);
				ypt = this.ptCloud{i}.Location(:, 2);

				% find closest pixel centre of each point
				x_idx = interp1(this.xc, 1:this.nx, xpt, 'nearest', 'extrap');
                y_idx = interp1(this.yc, 1:this.ny, ypt, 'nearest', 'extrap');
                this.ptCloud_pixels{i} = sub2ind([this.ny, this.nx], y_idx, x_idx); % pixel indices of all points in ptCloud

			end % for cycle end

			% compute metrics
			for i = 1:this.nPixels % iterate over all pixels of the mesh
				Z = [];
				groundPoints = 0;
				allPixelPoints = 0;

				fprintf("pixel: %d/%d\n", i, this.nPixels);

				% check whether i-th pixel should be computed
				if (~this.isPixelInCurve(i) && ~computeAllPixels)
					% if pixel is outside the curve and should not be computed -> continue to next pixel
					continue;
				end

				% find all points belonging to pixel "i"
				for j = 1:this.lasCount % iterate over all point clouds
					% find vegetation points
					isVegetation = this.ptAttributes{j}.Classification == this.VEG_LOW | ...
								   this.ptAttributes{j}.Classification == this.VEG_MED | ...
								   this.ptAttributes{j}.Classification == this.VEG_HIGH;
					% find points that belong to pixel "i"
					isInPixel = this.ptCloud_pixels{j} == i;

					% select only vegetation points belonging to pixel "i"
					selectedPoints = isInPixel & isVegetation;

					Z = [Z; this.ptCloud{j}.Location(selectedPoints, 3)];
					
					% find ground points which belong to pixel "i"
					% TODO: ci nestaci dat iba ~isVegetation, kedze body su bud iba vegetacia alebo ground
					isGround = this.ptAttributes{j}.Classification == this.GROUND;
					selectedPoints = isInPixel & isGround;
					
% 					if nnz(selectedPoints) ~= 0
% 						fprintf("...\n");
% 					end

					% find how many ground points there are within pixel "i"
					groundPoints = groundPoints + nnz(selectedPoints);
					% find how many points in total there are within pixel "i"
					allPixelPoints = allPixelPoints + nnz(isInPixel);
				end
				
				% if no points are found within pixel "i" -> continue to next pixel
				if allPixelPoints == 0
					continue;
				end
				
				% TODO: mozno zbytocna podmienka, asi staci len ta nad nou
				if isempty(Z)
					Z = 0;
				end

				% compute metrics
				maxZ     = max(Z);
				meanZ    = mean(Z);
				medianZ  = median(Z);
				p25Z     = prctile(Z, 25);
				p75Z     = prctile(Z, 75);
				p95Z     = prctile(Z, 95);

				PPR         = groundPoints / allPixelPoints;
				DAM_z       = nnz(Z > meanZ);
				BR_bellow_1 = nnz(Z < 1) / length(Z);
				BR_1_2      = nnz (Z > 1 & Z < 2) / length(Z);
				BR_2_3      = nnz (Z > 2 & Z < 3) / length(Z);
				BR_above_3  = nnz (Z > 3) / length(Z);
				BR_3_4      = nnz (Z > 3 & Z < 4) / length(Z);
				BR_4_5      = nnz (Z > 4 & Z < 5) / length(Z);
				BR_bellow_5 = nnz (Z < 5) / length(Z);
				BR_5_20     = nnz (Z > 5 & Z < 20) / length(Z);
				BR_above_20 = nnz (Z > 20) / length(Z);

				kurtZ    = kurtosis(Z);
				skewZ    = skewness(Z);

% 				if kurtZ > 7.5 || kurtZ < 1.5
% 
% 					fprintf("aaaa");
% 				end
				varZ     = var(Z);
				stdZ     = std(Z);
				coefVarZ = stdZ / meanZ;

				% save metrics to appropriate rasters
				this.metricsRasters.Hmax(i)     = maxZ;
				this.metricsRasters.Hmean(i)    = meanZ;
				this.metricsRasters.Hmedian(i)  = medianZ;
				this.metricsRasters.Hp25(i)     = p25Z;
				this.metricsRasters.Hp75(i)     = p75Z;
				this.metricsRasters.Hp95(i)     = p95Z;
				
				this.metricsRasters.PPR(i)		   = PPR;
				this.metricsRasters.DAM_z(i)	   = DAM_z;
				this.metricsRasters.BR_bellow_1(i) = BR_bellow_1;
				this.metricsRasters.BR_1_2(i)      = BR_1_2;
				this.metricsRasters.BR_2_3(i)      = BR_2_3;
				this.metricsRasters.BR_above_3(i)  = BR_above_3;
				this.metricsRasters.BR_3_4(i)      = BR_3_4;
				this.metricsRasters.BR_4_5(i)      = BR_4_5;
				this.metricsRasters.BR_bellow_5(i) = BR_bellow_5;
				this.metricsRasters.BR_5_20(i)     = BR_5_20;
				this.metricsRasters.BR_above_20(i) = BR_above_20;

				this.metricsRasters.Coeff_var_z(i) = coefVarZ;
				this.metricsRasters.Hkurt(i)       = kurtZ;
				this.metricsRasters.Hskew(i)       = skewZ;
				this.metricsRasters.Hstd(i)        = stdZ;
				this.metricsRasters.Hvar(i)        = varZ;

			end

			time = toc(time);
			fprintf("Compute rasters done in %0.3f s\n", time);

			temp = flipud(this.metricsRasters.Hmax);
			this.alphaData = ones(size(temp)); % vytvorenie pola, kde budu ulozene udaje o alpha pre kazdy pixel
			this.alphaData(isnan(temp)) = 0; % tam, kde su NaN hodnoty, tak budu priehladne

			timePassed = time;
		end % end of function computeMetricRasters

		function [this, timePassed] = computeRepresentativeMetrics(this)
			arguments
				this dataFeatureExtractorCurve
			end

			rasterNames = fieldnames(this.metricsRasters);
			reprMetrics = zeros(4, length(rasterNames));
			currentRaster = [];
			values = [];

			tic
			for i = 1:length(rasterNames)
				currentRaster = this.metricsRasters.(rasterNames{i});

				values = currentRaster(this.isPixelInCurve);

				reprMetrics(1,i) = mean(values);
				reprMetrics(2,i) = std(values);
				reprMetrics(3,i) = min(values);
				reprMetrics(4,i) = max(values);
			end

			this.representativeMetrics = (reprMetrics(:))';
			timePassed = toc;
		end

		function this = clipMetricRaster(this, options)
			arguments
				this 
				options.clipData (1,:) string {mustBeMember(options.clipData,{'Hmax', 'Hmean',...
					'Hmedian','Hp25','Hp75', 'Hp95'...
					'PPR','DAM_z','BR_bellow_1','BR_1_2','BR_2_3','BR_above_3','BR_3_4','BR_4_5','BR_bellow_5','BR_5_20','BR_above_20' ...
					'Coeff_var_z', 'Hkurt', 'Hskew', 'Hstd', 'Hvar'})}
				options.percentile (1,1) {mustBeNumeric, mustBeInRange(options.percentile, 0, 100)} = 97.5
			end
			
			chosenRaster = this.metricsRasters.(options.clipData);
			chosenPercentile = prctile(chosenRaster(:), options.percentile);

			chosenRaster(chosenRaster > chosenPercentile) = chosenPercentile;

			this.metricsRasters.(options.clipData) = chosenRaster;

			fprintf("Data " + options.clipData + " cliped\n");
		end

		function plotMetricRaster(this, options)
			arguments
				this dataFeatureExtractorCurve
				options.plotData (1,:) string {mustBeMember(options.plotData,["Hmax", "Hmean",...
					"Hmedian","Hp25","Hp75", "Hp95"...
					"PPR","DAM_z","BR_bellow_1","BR_1_2","BR_2_3","BR_above_3","BR_3_4","BR_4_5","BR_bellow_5","BR_5_20","BR_above_20" ...
					"Coeff_var_z", "Hkurt", "Hskew", "Hstd", "Hvar"])}
				options.clipPercentile (1,1) {mustBeNumeric, mustBeInRange(options.clipPercentile, 0, 100)} = 97.5
				options.plotCurve logical = 0
				options.plotMesh logical = 0
			end

			% save copy of raster to be ploted
			chosenRaster = this.metricsRasters.(options.plotData);
			chosenRaster = flipud(chosenRaster);

			if ismember(options.plotData, ["Coeff_var_z", "Hkurt", "Hskew", "Hvar"])
				% compute chosen percentile for clipping
				chosenPercentile = prctile(chosenRaster(:), options.clipPercentile);

				% clip raster data according to the chosen percentile
				chosenRaster(chosenRaster > chosenPercentile) = chosenPercentile;
			end
			
			imagesc([this.xc(1) this.xc(end)], [this.yc(end) this.yc(1)], chosenRaster,...
				'AlphaData', this.alphaData)
			title(this.plotTitles.(options.plotData))

			hold on
			if (options.plotMesh)
% 				plot(this.curve, "EdgeColor", "magenta", "FaceAlpha", 0,"LineWidth",3,"LineStyle","-")
				this.plotMesh();
			end
			hold off
			axis equal
			axis xy

		end % end of function plotMetricRaster

		function plotAllMetricRasters(this, options)
			arguments
				this
				options.colormap (1,:) string {mustBeMember(options.colormap,{'parula','turbo','hsv',...
					'hot','cool','spring','summer','autumn','winter','gray','bone','copper',...
					'pink','jet','lines','colorcube','prism','flag','white'})} = 'jet'
				options.plotCurves = {}
			end
 
			rasterNames = fieldnames(this.metricsRasters);
			
			for i = numel(rasterNames):-1:1
				figure('WindowState','maximized')
				plotMetricRaster(this,"plotData",rasterNames{i});
				colormap(options.colormap)
				colorbar
				axis equal
				axis xy

				if ~isempty(options.plotCurves)
					plotCurves(options.plotCurves,"Color",'magenta','MarkerSize',20)
				end
			end

		end % end of function plotAllMetricRasters

		function exportMetricRaster(this, fileName, options)
			arguments
				this 
				fileName (1,:) string
				options.exportLayer (1,:) string {mustBeMember(options.exportLayer,{'Hmax', 'Hmean',...
					'Hmedian','Hp25','Hp75', 'Hp95'...
					'PPR','DAM_z','BR_bellow_1','BR_1_2','BR_2_3','BR_above_3','BR_3_4','BR_4_5','BR_bellow_5','BR_5_20','BR_above_20' ...
					'Coeff_var_z', 'Hkurt', 'Hskew', 'Hstd', 'Hvar'})}
			end

			exportName = this.exportTitles.(options.exportLayer) + fileName;
			geotiffwrite(exportName, this.metricsRasters.(options.exportLayer),...
				this.RR,"CoordRefSysCode",8353);
			
			fprintf("Metric " + options.exportLayer + " exported\n");
		end

		function exportAllMetricRasters(this, fileName)
			arguments
				this 
				fileName (1,:) string
			end

			rasterNames = fieldnames(this.metricsRasters);
			
			for i = 1:numel(rasterNames)
				exportName = this.exportTitles.(rasterNames{i}) + fileName;
				
				geotiffwrite(exportName, this.metricsRasters.(rasterNames{i}),...
				this.RR,"CoordRefSysCode",8353);
			end
		
			fprintf("All metrics exported\n");
		end
		
	end
end

