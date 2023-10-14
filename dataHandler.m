classdef dataHandler
	% popis triedy
	properties
		% DTM stuff
		DTM								% DTM data
		RasterReference					% RasterReference for DTM data
		alphaData_DTM					% set to zero for pixels with NaN values

		% 2D mesh parameters
		x1, x2, y1, y2                  % x/y pixel centers range
        nx, ny                          % number of pixels in x/y direction
        hx, hy                          % pixel sizes
        x, y                            % pixel edge coordinates
        xc, yc, Xc, Yc                  % pixel center coordinates (and meshgrids)
        
		% PC stuff
		ptCloud							% cell of all point clouds
		ptCloud_normalized				% cell of normalized point cloud (only vegetation)
		lasCount						% number of point clouds
		ptAttributes					% cell of classification for all points clouds
		ptCloud_pixels					% DTM pixel indices of all points in ptCloud

		metricsRasters
		
		GROUND = 2
		VEG_LOW = 3
		VEG_MED = 4
		VEG_HIGH = 5
	
	end % end of properties

	methods
		% CONSTRUCTOR
		function this = dataHandler(DTM, RR, ptCloud, ptAttributes)
% 			arguments
% 				H (:,:) double
% 				R (1,1) map.rasterref.MapCellsReference
% 				ptCloud 
% 				ptAttributes
% 			end

			this.DTM = DTM;
			this.RasterReference = RR;
			this.ptCloud = ptCloud;
			this.ptAttributes = ptAttributes;

% 			this.hx = round(RR.CellExtentInWorldX);
% 			this.hy = round(RR.CellExtentInWorldY);
% 
% 			this.nx = round(RR.RasterExtentInWorldX);
%             this.ny = round(RR.RasterExtentInWorldY);

			this.hx = RR.CellExtentInWorldX;
			this.hy = RR.CellExtentInWorldY;

			this.nx = RR.RasterExtentInWorldX;
            this.ny = RR.RasterExtentInWorldY;

			this.x1 = RR.XWorldLimits(1);
            this.x2 = this.x1 + (this.nx-1)*this.hx;
            this.y1 = RR.YWorldLimits(1);
            this.y2 = this.y1 + (this.ny-1)*this.hy;

			this.alphaData_DTM = ones(size(this.DTM)); % vytvorenie pola, kde budu ulozene udaje o alpha pre kazdy pixel
			this.alphaData_DTM(isnan(this.DTM)) = 0; % tam, kde su NaN hodnoty, tak budu priehladne
			
			this.metricsRasters = struct;

		end % end of constructor
		
		% MESH
		function this = meshPlane(this, padding)
			arguments
				this (1,1) dataHandler
				padding	(1,1) double = 0
			end

			if padding > 0
                this.DTM = padarray(this.DTM, [padding,padding], NaN, 'both');
                
                this.nx = this.nx + 2 * padding;
                this.ny = this.ny + 2 * padding;
                this.x1 = this.x1 - padding * this.hx;
                this.x2 = this.x2 + padding * this.hx;
                this.y1 = this.y1 - padding * this.hy;
                this.y2 = this.y2 + padding * this.hy;
			end

			this.x = this.x1 - this.hx/2:this.hx:this.x2 + this.hx/2;
            this.y = this.y1 - this.hy/2:this.hy:this.y2 + this.hy/2;
            
            % specification of pixel centers (interpolation points for interp1)
            this.xc = linspace(this.x1, this.x2, this.nx);
            this.yc = linspace(this.y1, this.y2, this.ny);
            [this.Xc, this.Yc] = meshgrid(this.xc, this.yc);

			fprintf("Mesh plane done\n");

		end % end of mesh plane

		function plotMesh(this, options)
			arguments
				this (1,1) dataHandler
				options.Color = [0.3, 0.3, 0.3]
				options.LineWidth (1,1) {mustBeNumeric} = 4
			end
			
			hold on
			for i = 1:length(this.x)
				plot([this.x(i) this.x(i)],[this.y(1) this.y(end)],...
					'Color', options.Color, 'LineWidth', options.LineWidth) %y grid lines
			end

			for i = 1:length(this.y)
				plot([this.x(1) this.x(end)],[this.y(i) this.y(i)],...
					'Color', options.Color, 'LineWidth', options.LineWidth) %x grid lines
			end
			hold off
			
		end % end of function plotMesh

		function this = computeTerrainAttributes(this)
			% najdenie pixelov, kde su NaN hodnoty, aby sa nevykreslovali 
			this.alphaData_DTM = ones(size(this.DTM)); % vytvorenie pola, kde budu ulozene udaje o alpha pre kazdy pixel
			this.alphaData_DTM(isnan(this.DTM)) = 0; % tam, kde su NaN hodnoty, tak budu priehladne

		end % end of function computeTerrainAttributes

		function plotTerrain(this, options)
			arguments
				this (1,1) dataHandler
				options.ncolors double = 256
			end

			imagesc([this.x1 this.x2], [this.y1 this.y2], this.DTM,...
					'AlphaData', this.alphaData_DTM)
			demcmap(this.DTM, options.ncolors);
			axis equal
			axis xy

		end % end of function plotTerrain

		function plotTerrain3D(this, options)
			arguments
				this (1,1) dataHandler
				options.useData (1,:) string {mustBeMember(options.useData,{'H', 'H_filt', 'H_subtracted'})} = 'H_filt'
				options.ncolors double {mustBePositive} = 256
				options.EdgeColor = 'k'
				options.EdgeAlpha {mustBeNumeric, mustBeGreaterThanOrEqual(options.EdgeAlpha, 0),...
					mustBeLessThanOrEqual(options.EdgeAlpha, 1)} = 0.05 % musi byt 0 <= alpha <= 1
				options.FaceAlpha {mustBeNumeric, mustBeGreaterThanOrEqual(options.FaceAlpha, 0),...
					mustBeLessThanOrEqual(options.FaceAlpha, 1)} = 1.0 % musi byt 0 <= alpha <= 1
			end

			s = surf(this.Xc, this.Yc, this.DTM);
			demcmap(this.DTM, options.ncolors);

			s.EdgeColor = options.EdgeColor;
			s.EdgeAlpha = options.EdgeAlpha;
			s.FaceAlpha = options.FaceAlpha;
			
		end % end of function plotTerrain3D

		function this = computePointCloudAttributes(this)
			% kontrola, ci je priradeny nejaky point cloud
			if isempty(this.ptCloud)
				warning('No point cloud available');
				return;
			end

			fprintf('Computing plane point cloud attributes... ');
            time = tic;

			this.lasCount = length(this.ptCloud);
			this.ptCloud_pixels = cell(this.lasCount, 1);
% 			this.pixel_ptCloud_indices = cell(this.ny, this.nx); % empty cell array

			for i = 1:this.lasCount
				xpt = this.ptCloud{i}.Location(:, 1);
				ypt = this.ptCloud{i}.Location(:, 2);

				% find closest pixel centre of each point
				x_idx = interp1(this.xc, 1:this.nx, xpt, 'nearest', 'extrap');
                y_idx = interp1(this.yc, 1:this.ny, ypt, 'nearest', 'extrap');
                this.ptCloud_pixels{i} = sub2ind([this.ny, this.nx], y_idx, x_idx); % pixel indices of all points in ptCloud

				% find indices of all ptCloud points in each pixel of the mesh
%                 for j = 1:this.ptCloud{i}.Count % for all points in ptCloud
%                     current = this.pixel_ptCloud_indices{this.ptCloud_pixels{i}(j)};
%                     this.pixel_ptCloud_indices{this.ptCloud_pixels{i}(j)} = [current; [i,j]]; % add point j to pixel this.ptCloud_pixels(i,:)
%                 end

			end % for cycle end

			% filter points belonging to NaN pixels of DTM
			for i = 1:this.lasCount
				selectedPoints = ~isnan(this.DTM(this.ptCloud_pixels{i}));
				
				% select only points from non-NaN pixels, attributes and pixel indices
				this.ptCloud{i} = select(this.ptCloud{i}, selectedPoints);
				this.ptAttributes{i} = lidarPointAttributes("Classification",...
									    this.ptAttributes{i}.Classification(selectedPoints));
				this.ptCloud_pixels{i} = this.ptCloud_pixels{i}(selectedPoints);

			end

			time = toc(time);
			fprintf('done in %0.2f s\n', time)

		end % end of function computePointCloudAttributes

		function this = normalizePtCloud(this, options)
			arguments
				this (1,1) dataHandler
				options.method (1,:) string {mustBeMember(options.method,{'DTM', 'LowestPoint'})} = 'DTM'
			end
			
			time = tic;
			this.ptCloud_normalized = cell(this.lasCount, 1);
			
			if strcmp(options.method, "DTM")
				for i = 1:this.lasCount
% 					vegetation = ismember(this.ptAttributes{i}.Classification, [this.VEG_LOW, this.VEG_MED, this.VEG_HIGH]);
					
					X = this.ptCloud{i}.Location(:, 1);
					Y = this.ptCloud{i}.Location(:, 2);
					Z = this.ptCloud{i}.Location(:, 3);
					
					Z = Z - this.DTM(this.ptCloud_pixels{i});

					Z(Z < 0) = 0;
	
					this.ptCloud_normalized{i} = pointCloud([X, Y, Z]);
				end
			end

			if strcmp(options.method, "LowestPoint")
				
				for i = 1:this.lasCount % iterate over all point clouds
% 					vegetation = ismember(this.ptAttributes{i}.Classification, [this.VEG_LOW, this.VEG_MED, this.VEG_HIGH]);
					
					X = this.ptCloud{i}.Location(:, 1);
					Y = this.ptCloud{i}.Location(:, 2);
					Z = this.ptCloud{i}.Location(:, 3);
					% TODO: pridat "Intensity" do normalizovaneho point cloudu

					for j = 1:length(this.DTM(:)) % iterate through each pixel of mesh
						pixelPoints = this.ptCloud_pixels{i}(:) == j;
						if nnz(pixelPoints) == 0
							continue;
						end
						
						Z(pixelPoints) = Z(pixelPoints) - min(Z(pixelPoints));
					end
		
					this.ptCloud_normalized{i} = pointCloud([X, Y, Z]);
				end
			end
			
			time = toc(time);
			
			fprintf("Normalize points done in %0.3f s\n", time);
		end % end of normalize point cloud
		
		function this = computeMetricRasters(this)
			% allocate space for metric rasters
			% ECOSYSTEM HEIGHT
			this.metricsRasters.Hmax = NaN(this.ny, this.nx);
			this.metricsRasters.Hmean = NaN(this.ny, this.nx);
			this.metricsRasters.Hmedian = NaN(this.ny, this.nx);
			this.metricsRasters.Hp25 = NaN(this.ny, this.nx);
			this.metricsRasters.Hp75 = NaN(this.ny, this.nx);
			this.metricsRasters.Hp95 = NaN(this.ny, this.nx);
			NaN
			% ECOSYSTEM COVER
			this.metricsRasters.PPR = NaN(this.ny, this.nx);

			% ECOSYSTEM STRUCTURAL COMPLEXITY
			this.metricsRasters.Coeff_var_z = NaN(this.ny, this.nx);
			this.metricsRasters.Hkurt = NaN(this.ny, this.nx);
			this.metricsRasters.Hskew = NaN(this.ny, this.nx);
			this.metricsRasters.Hstd = NaN(this.ny, this.nx);
			this.metricsRasters.Hvar = NaN(this.ny, this.nx);
				
			time = tic;

			for i = 1:length(this.DTM(:)) % iterate over all pixels
				% skip pixels of DTM where there is NaN value
				if isnan(this.DTM(i))
					continue;
				end
				
				Z = [];
				groundPoints = 0;
				allPixelPoints = 0;

				fprintf("pixel: %d/%d\n", i, length(this.DTM(:)))
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

					Z = [Z; this.ptCloud_normalized{j}.Location(selectedPoints, 3)];
					
					% find ground points which belong to pixel "i"
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

				PPR      = groundPoints / allPixelPoints;

% 				if (groundPoints > allPixelPoints)
% 					fprintf("aaaa\n");
% 				end

				kurtZ    = kurtosis(Z);
				skewZ    = skewness(Z);
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
				
				this.metricsRasters.PPR(i)		 = PPR;

				this.metricsRasters.Coeff_var_z(i) = coefVarZ;
				this.metricsRasters.Hkurt(i)       = kurtZ;
				this.metricsRasters.Hskew(i)       = skewZ;
				this.metricsRasters.Hstd(i)        = stdZ;
				this.metricsRasters.Hvar(i)        = varZ;

			end

			time = toc(time);
			fprintf("Compute rasters done in %0.3f s\n", time);
		end

		function plotMetricRaster(this, options)
			arguments
				this (1,1) dataHandler
				options.plotData (1,:) string {mustBeMember(options.plotData,{'Hmax', 'Hmean',...
					'Hmedian','Hp25','Hp75', 'Hp95',...
					'PPR', ...
					'Coeff_var_z', 'Hkurt', 'Hskew', 'Hstd', 'Hvar'})}
			end
			
			if strcmp(options.plotData, 'Hmax')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmax,...
					'AlphaData', this.alphaData_DTM)
				title 'Max of normalized height'
			end

			if strcmp(options.plotData, 'Hmean')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmean,...
					'AlphaData', this.alphaData_DTM)
				title 'Mean of normalized height'
			end

			if strcmp(options.plotData, 'Hmedian')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmedian,...
					'AlphaData', this.alphaData_DTM)
				title 'Median of normalized height'
			end

			if strcmp(options.plotData, 'Hp25')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp25,...
					'AlphaData', this.alphaData_DTM)
				title '25th percentile of normalized height'
			end

			if strcmp(options.plotData, 'Hp75')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp75,...
					'AlphaData', this.alphaData_DTM)
				title '75th percentile of normalized height'
			end

			if strcmp(options.plotData, 'Hp95')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp95,...
					'AlphaData', this.alphaData_DTM)
				title '95th percentile of normalized height'
			end

			if strcmp(options.plotData, 'PPR')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.PPR,...
					'AlphaData', this.alphaData_DTM)
				title 'Pulse penetration ratio'
			end

			if strcmp(options.plotData, 'Coeff_var_z')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Coeff_var_z,...
					'AlphaData', this.alphaData_DTM)
				title 'Coefficient of variation'
			end

			if strcmp(options.plotData, 'Hkurt')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hkurt,...
					'AlphaData', this.alphaData_DTM)
				title 'Kurtosis of normalized height'
			end

			if strcmp(options.plotData, 'Hskew')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hskew,...
					'AlphaData', this.alphaData_DTM)
				title 'Skewness of normalized height'
			end

			if strcmp(options.plotData, 'Hstd')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hstd,...
					'AlphaData', this.alphaData_DTM)
				title 'Standard deviation of normalized height'
			end

			if strcmp(options.plotData, 'Hvar')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hvar,...
					'AlphaData', this.alphaData_DTM)
				title 'Variance of normalized height'
			end

			colormap jet
			colorbar
			axis equal
			axis xy

		end % end of function plotMetricRaster

		function plotMetricRasters(this)
			arguments
				this (1,1) dataHandler
			end

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmax,...
				'AlphaData', this.alphaData_DTM)
			title 'Max of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmean,...
				'AlphaData', this.alphaData_DTM)
			title 'Mean of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmedian,...
				'AlphaData', this.alphaData_DTM)
			title 'Median of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy


			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp25,...
				'AlphaData', this.alphaData_DTM)
			title '25th percentile of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp75,...
				'AlphaData', this.alphaData_DTM)
			title '75th percentile of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp95,...
				'AlphaData', this.alphaData_DTM)
			title '95th percentile of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.PPR,...
				'AlphaData', this.alphaData_DTM)
			title 'Pulse penetration ratio'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Coeff_var_z,...
				'AlphaData', this.alphaData_DTM)
			title 'Coefficient of variation'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hkurt,...
				'AlphaData', this.alphaData_DTM)
			title 'Kurtosis of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hskew,...
				'AlphaData', this.alphaData_DTM)
			title 'Skewness of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hstd,...
				'AlphaData', this.alphaData_DTM)
			title 'Standard deviation of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

			figure
			imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hvar,...
				'AlphaData', this.alphaData_DTM)
			title 'Variance of normalized height'
			colormap jet
			colorbar
			axis equal
			axis xy

		end % end of function plotMetricRasters

		function exportMetricRaster(this, fileName, options)
			arguments
				this dataHandler
				fileName (1,:) string
				options.exportLayer (1,:) string {mustBeMember(options.exportLayer,{'Hmax', 'Hmean',...
					'Hmedian','Hp25','Hp75', 'Hp95',...
					'PPR', ...
					'Coeff_var_z', 'Hkurt', 'Hskew', 'Hstd', 'Hvar'})}
			end

			if strcmp(options.exportLayer, 'Hmax')
				geotiffwrite("1_Hmax_" + fileName, this.metricsRasters.Hmax,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hmean')
				geotiffwrite("2_Hmean_" + fileName, this.metricsRasters.Hmean,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hmedian')
				geotiffwrite("3_Hmedian_" + fileName, this.metricsRasters.Hmedian,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hp25')
				geotiffwrite("4_Hp25_" + fileName, this.metricsRasters.Hp25,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hp75')
				geotiffwrite("6_Hp75_" + fileName, this.metricsRasters.Hp75,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hp95')
				geotiffwrite("7_Hp95_" + fileName, this.metricsRasters.Hp95,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'PPR')
				geotiffwrite("8_PPR_" + fileName, this.metricsRasters.PPR,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Coeff_var_z')
				geotiffwrite("19_Coeff_var_z_" + fileName, this.metricsRasters.Coeff_var_z,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hkurt')
				geotiffwrite("21_Hkurt_" + fileName, this.metricsRasters.Hkurt,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hskew')
				geotiffwrite("23_Hskew_" + fileName, this.metricsRasters.Hskew,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hstd')
				geotiffwrite("24_Hstd_" + fileName, this.metricsRasters.Hstd,...
					this.RasterReference,"CoordRefSysCode",8353);
			end

			if strcmp(options.exportLayer, 'Hvar')
				geotiffwrite("25_Hvar_" + fileName, this.metricsRasters.Hvar,...
					this.RasterReference,"CoordRefSysCode",8353);
			end
			
			fprintf("Metric " + options.exportLayer + " exported\n");
		end

		function exportAllMetricRasters(this, fileName)
			arguments
				this dataHandler
				fileName (1,:) string
			end
			
			% Maximum of normalized height
			geotiffwrite("1_Hmax_" + fileName, this.metricsRasters.Hmax,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Mean of normalized height
			geotiffwrite("2_Hmean_" + fileName, this.metricsRasters.Hmean,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Median of normalized height
			geotiffwrite("3_Hmedian_" + fileName, this.metricsRasters.Hmedian,...
				this.RasterReference,"CoordRefSysCode",8353);
			% 25th percentile of normalized height
			geotiffwrite("4_Hp25_" + fileName, this.metricsRasters.Hp25,...
				this.RasterReference,"CoordRefSysCode",8353);
			% 75th percentile of normalized height
			geotiffwrite("6_Hp75_" + fileName, this.metricsRasters.Hp75,...
				this.RasterReference,"CoordRefSysCode",8353);
			% 95th percentile of normalized height
			geotiffwrite("7_Hp95_" + fileName, this.metricsRasters.Hp95,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Pulse Penetration Ratio
			geotiffwrite("8_PPR_" + fileName, this.metricsRasters.PPR,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Coefficient of variation
			geotiffwrite("19_Coeff_var_z_" + fileName, this.metricsRasters.Coeff_var_z,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Kurtosis of normalized height
			geotiffwrite("21_Hkurt_" + fileName, this.metricsRasters.Hkurt,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Skewness of normalized height
			geotiffwrite("23_Hskew_" + fileName, this.metricsRasters.Hskew,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Standard deviation of normalized height
			geotiffwrite("24_Hstd_" + fileName, this.metricsRasters.Hstd,...
				this.RasterReference,"CoordRefSysCode",8353);
			% Variance of normalized height
			geotiffwrite("25_Hvar_" + fileName, this.metricsRasters.Hvar,...
				this.RasterReference,"CoordRefSysCode",8353);

			fprintf("All metrics exported\n");
		end

		function plotPtCloud2D(this, colorData, selectedClass) 
			arguments
				this dataHandler
				colorData 
				selectedClass 
			end
			
			
			if this.lasCount == 0
				warning('No Point Cloud')
				return
			end
			
% 			figure
			% TH: pri vykreslovani vsetkych bodov nie je funkcia scatter zrovna najlepsia, co sa tyka plynulosti 
% 			for i = 1:lasCount
% 				scatter(this.ptCloud{i}.Location(:, 1), this.ptCloud{i}.Location(:, 2), 1, colorData{i}/255)
% 				hold on
% 			end
			% TH: preto sa vykresluje zatial vzdy iba jedna vybrana trieda
			for i = 1:this.lasCount
				classMember = ismember(this.ptAttributes{i}.Classification, selectedClass);
				if any(classMember)
					scatter(this.ptCloud{i}.Location(classMember, 1), this.ptCloud{i}.Location(classMember, 2),...
						10, colorData{i}(classMember, :)/255, 'filled')
					%hold on
				else
					%warning('Nothing to plot')
				end
			end
			
			%hold off
            %title 'Terrain Point Cloud'
%             axis ij
            %axis equal
		end % end of plot point cloud 2D

		function plotPtCloud3D(this, colorMap, selectedClasses, options)
            arguments
				this (1,1) dataHandler
				colorMap (:, 3)
				selectedClasses
				options.useData (1,:) string {mustBeMember(options.useData,{'original', 'normalized'})} = 'original'
			end

			if this.lasCount == 0
				warning('Nothing to plot')
				return
			end
        
			% plot only points from selected classes
			hold on

			if strcmp(options.useData, 'original')

				for i = 1:this.lasCount
					classMember = ismember(this.ptAttributes{i}.Classification, selectedClasses);
					if any(classMember)
						colorData = reshape(label2rgb(this.ptAttributes{i}.Classification, colorMap, 'k'), [], 3);
						pcshow(this.ptCloud{i}.Location(classMember, :), colorData(classMember, :))
					else
						warning('Nothing to plot')
					end
				end
			
			elseif strcmp(options.useData, 'normalized')
				
				for i = 1:this.lasCount
					classMember = ismember(this.ptAttributes{i}.Classification, selectedClasses);
					if any(classMember)
						colorData = reshape(label2rgb(this.ptAttributes{i}.Classification, colorMap, 'k'), [], 3);
						pcshow(this.ptCloud_normalized{i}.Location(classMember, :), colorData(classMember, :))
					else
						warning('Nothing to plot')
					end
				end


			end
			
			hold off
			
		end % end of plot point cloud 3D

	end % end of methods
end % end of classdef