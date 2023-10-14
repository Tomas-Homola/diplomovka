classdef dataFeatureExtractor
	% Class for feature extraction from the given point cloud
	% 1) load normalized point cloud from dataPreprocessor and RasterReference
	% 2) create MESH with pixel size h x h meters
	% 3) find indices within the new mesh
	% 4) compute metrics and save into preallocated space
	
	properties
		RR_orig map.rasterref.MapCellsReference % RasterReference of original DTM
		RR_new map.rasterref.MapCellsReference % RasterReference of the new mesh

		% 2D mesh parameters
		x1, x2, y1, y2                  % x/y pixel centers range
        nx, ny                          % number of pixels in x/y direction
        h                               % desired pixel size
        x, y                            % pixel edge coordinates
        xc, yc, Xc, Yc                  % pixel center coordinates (and meshgrids)
		numOfMeshPixels
		alphaData

		% Point cloud stuff
		lasCount						% number of point clouds
		ptCloud_norm				    % cell of normalized point cloud
		ptAttributes_norm				% cell of classification for all points clouds
		ptCloud_norm_pixels				% pixel indices of all points within the mesh

		% Important point cloud classes
		GROUND = 2
		VEG_LOW = 3
		VEG_MED = 4
		VEG_HIGH = 5

		% structure for rasters of individual metrics computed from point cloud
		metricsRasters
	end
	
	methods
		function this = dataFeatureExtractor(RR_orig, ptCloud_norm, ptAttributes_norm, h)
			arguments
				RR_orig map.rasterref.MapCellsReference
				ptCloud_norm
				ptAttributes_norm
				h (1,1) {mustBeNumeric, mustBeInteger, mustBeGreaterThanOrEqual(h, 1)}
			end

			this.RR_orig = RR_orig;
			this.ptCloud_norm = ptCloud_norm;
			this.ptAttributes_norm = ptAttributes_norm;
			this.lasCount = length(ptCloud_norm);

			this.h = h;

			nxOld = RR_orig.RasterExtentInWorldX;
            nyOld = RR_orig.RasterExtentInWorldY;

			% compute new raster extent based on the desired pixel size
			this.nx = (nxOld - mod(nxOld, h)) / h;
			this.ny = (nyOld - mod(nyOld, h)) / h;

			this.numOfMeshPixels = this.ny * this.nx;
			
			this.x1 = RR_orig.XWorldLimits(1);
            this.x2 = this.x1 + (this.nx-1)*this.h; % ci to je dobre?
            this.y1 = RR_orig.YWorldLimits(1);
            this.y2 = this.y1 + (this.ny-1)*this.h;

			% create new RasterReference object for the new mesh

		end % end of constructor

		function this = meshPlane(this, padding)
			arguments
				this (1,1) dataFeatureExtractor
				padding	(1,1) double = 0
			end

% 			if padding > 0
%                 this.DTM = padarray(this.DTM, [padding,padding], NaN, 'both');
%                 
%                 this.nx = this.nx + 2 * padding;
%                 this.ny = this.ny + 2 * padding;
%                 this.x1 = this.x1 - padding * this.hx;
%                 this.x2 = this.x2 + padding * this.hx;
%                 this.y1 = this.y1 - padding * this.hy;
%                 this.y2 = this.y2 + padding * this.hy;
% 			end
	
			this.x = this.x1 - this.h/2:this.h:this.x2 + this.h/2;
            this.y = this.y1 - this.h/2:this.h:this.y2 + this.h/2;
            
            % specification of pixel centers (interpolation points for interp1)
            this.xc = linspace(this.x1, this.x2, this.nx);
            this.yc = linspace(this.y1, this.y2, this.ny);
            [this.Xc, this.Yc] = meshgrid(this.xc, this.yc);

			fprintf("Mesh plane done\n");

		end % end of mesh plane

		function plotMesh(this, options)
			arguments
				this (1,1) dataFeatureExtractor
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

			axis equal
			
		end % end of function plotMesh

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

			% find pixel indices for all points within mesh
			this.ptCloud_norm_pixels = cell(this.lasCount, 1);
			for i = 1:this.lasCount
				xpt = this.ptCloud_norm{i}.Location(:, 1);
				ypt = this.ptCloud_norm{i}.Location(:, 2);

				% find closest pixel centre of each point
				x_idx = interp1(this.xc, 1:this.nx, xpt, 'nearest', 'extrap');
                y_idx = interp1(this.yc, 1:this.ny, ypt, 'nearest', 'extrap');
                this.ptCloud_norm_pixels{i} = sub2ind([this.ny, this.nx], y_idx, x_idx); % pixel indices of all points in ptCloud

			end % for cycle end

			% compute metrics
			for i = 1:this.numOfMeshPixels % iterate over all pixels of the mesh
				Z = [];
				groundPoints = 0;
				allPixelPoints = 0;

				fprintf("pixel: %d/%d\n", i, this.numOfMeshPixels)
				% find all points belonging to pixel "i"
				for j = 1:this.lasCount % iterate over all point clouds
					% find vegetation points
					isVegetation = this.ptAttributes_norm{j}.Classification == this.VEG_LOW | ...
								   this.ptAttributes_norm{j}.Classification == this.VEG_MED | ...
								   this.ptAttributes_norm{j}.Classification == this.VEG_HIGH;
					% find points that belong to pixel "i"
					isInPixel = this.ptCloud_norm_pixels{j} == i;

					% select only vegetation points belonging to pixel "i"
					selectedPoints = isInPixel & isVegetation;

					Z = [Z; this.ptCloud_norm{j}.Location(selectedPoints, 3)];
					
					% find ground points which belong to pixel "i"
					isGround = this.ptAttributes_norm{j}.Classification == this.GROUND;
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

			this.alphaData = ones(size(this.metricsRasters.Hmax)); % vytvorenie pola, kde budu ulozene udaje o alpha pre kazdy pixel
			this.alphaData(isnan(this.metricsRasters.Hmax)) = 0; % tam, kde su NaN hodnoty, tak budu priehladne
		end % end of function computeMetricRasters

		function plotMetricRaster(this, options)
			arguments
				this (1,1) dataFeatureExtractor
				options.plotData (1,:) string {mustBeMember(options.plotData,{'Hmax', 'Hmean',...
					'Hmedian','Hp25','Hp75', 'Hp95',...
					'PPR', ...
					'Coeff_var_z', 'Hkurt', 'Hskew', 'Hstd', 'Hvar'})}
			end
			
			if strcmp(options.plotData, 'Hmax')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmax,...
					'AlphaData', this.alphaData)
				title 'Max of normalized height'
			end

			if strcmp(options.plotData, 'Hmean')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmean,...
					'AlphaData', this.alphaData)
				title 'Mean of normalized height'
			end

			if strcmp(options.plotData, 'Hmedian')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hmedian,...
					'AlphaData', this.alphaData)
				title 'Median of normalized height'
			end

			if strcmp(options.plotData, 'Hp25')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp25,...
					'AlphaData', this.alphaData)
				title '25th percentile of normalized height'
			end

			if strcmp(options.plotData, 'Hp75')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp75,...
					'AlphaData', this.alphaData)
				title '75th percentile of normalized height'
			end

			if strcmp(options.plotData, 'Hp95')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hp95,...
					'AlphaData', this.alphaData)
				title '95th percentile of normalized height'
			end

			if strcmp(options.plotData, 'PPR')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.PPR,...
					'AlphaData', this.alphaData)
				title 'Pulse penetration ratio'
			end

			if strcmp(options.plotData, 'Coeff_var_z')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Coeff_var_z,...
					'AlphaData', this.alphaData)
				title 'Coefficient of variation'
			end

			if strcmp(options.plotData, 'Hkurt')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hkurt,...
					'AlphaData', this.alphaData)
				title 'Kurtosis of normalized height'
			end

			if strcmp(options.plotData, 'Hskew')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hskew,...
					'AlphaData', this.alphaData)
				title 'Skewness of normalized height'
			end

			if strcmp(options.plotData, 'Hstd')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hstd,...
					'AlphaData', this.alphaData)
				title 'Standard deviation of normalized height'
			end

			if strcmp(options.plotData, 'Hvar')
				imagesc([this.x1 this.x2], [this.y1 this.y2], this.metricsRasters.Hvar,...
					'AlphaData', this.alphaData)
				title 'Variance of normalized height'
			end

			colormap jet
			colorbar
			axis equal
			axis xy

		end % end of function plotMetricRaster
		
	end
end

