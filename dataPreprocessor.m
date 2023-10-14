classdef dataPreprocessor
	% Class for preprocessing point cloud data:
	% 1) load DTM and original point cloud ptCloud_0
	% 2) creates mesh identical to DTM, "meshPlane()"
	% 3) filter point cloud ptCloud_0:
	%   3.1) select points within DTM bounds
	%	3.2) find indices of points within mesh
	%   3.3) select only vegetation and ground points
	%   3.4) remove points within NaN pixels of DTM (necessary if normalization via DTM)
	% 4) normalize point cloud PC_0 using DTM or min(Z) within pixel -> PC_n

	properties
		% DTM stuff
		DTM								   % DTM data
		RR map.rasterref.MapCellsReference % RasterReference for DTM data

		% 2D mesh parameters
		x1, x2, y1, y2                  % x/y pixel centers range
        nx, ny                          % number of pixels in x/y direction
        hx, hy                          % pixel sizes
        x, y                            % pixel edge coordinates
        xc, yc, Xc, Yc                  % pixel center coordinates (and meshgrids)
        
		% PC stuff
		ptCloud_0						% cell of all point clouds
		ptCloud_norm				    % cell of normalized point cloud
		lasCount						% number of point clouds
		ptAttributes					% cell of classification for all points clouds
		ptCloud_0_pixels				% DTM pixel indices of all points in ptCloud
		ptAttributes_norm
		ptCloud_norm_pixels

		GROUND = 2
		VEG_LOW = 3
		VEG_MED = 4
		VEG_HIGH = 5
	
	end % end of properties

	methods
		% CONSTRUCTOR
		function this = dataPreprocessor(DTM, RR, ptCloud, ptAttributes)

			this.DTM = DTM;
			this.RR = RR;
			this.ptCloud_0 = ptCloud;
			this.ptAttributes = ptAttributes;

			this.hx = RR.CellExtentInWorldX;
			this.hy = RR.CellExtentInWorldY;

			this.nx = RR.RasterExtentInWorldX;
            this.ny = RR.RasterExtentInWorldY;

			this.x1 = RR.XWorldLimits(1);
            this.x2 = this.x1 + (this.nx-1)*this.hx;
            this.y1 = RR.YWorldLimits(1);
            this.y2 = this.y1 + (this.ny-1)*this.hy;

		end % end of constructor
		
		% MESH
		function this = meshPlane(this, padding)
			arguments
				this (1,1) dataPreprocessor
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

		function this = filterPointCloud(this)
			% kontrola, ci je priradeny nejaky point cloud
			if isempty(this.ptCloud_0)
				error('No point cloud available');
			end

			fprintf('Computing plane point cloud attributes... ');
            time = tic;

			this.lasCount = length(this.ptCloud_0);
			this.ptCloud_norm = cell(this.lasCount, 1);
			% 3.1) select points within DTM bounds
			for i = 1:this.lasCount
				selectedPoints = this.ptCloud_0{i}.Location(:,1) > this.RR.XWorldLimits(1) &...
					             this.ptCloud_0{i}.Location(:,1) < this.RR.XWorldLimits(2) &...
							     this.ptCloud_0{i}.Location(:,2) > this.RR.YWorldLimits(1) &...
								 this.ptCloud_0{i}.Location(:,2) < this.RR.XWorldLimits(2);

				% select only points within DTM bounds
				this.ptCloud_norm{i}      = select(this.ptCloud_0{i}, selectedPoints);
				this.ptAttributes_norm{i} = lidarPointAttributes("Classification",...
									        this.ptAttributes{i}.Classification(selectedPoints));
			end
			
			% 3.2) find indices of points within mesh
			this.ptCloud_norm_pixels = cell(this.lasCount, 1);
			for i = 1:this.lasCount
				xpt = this.ptCloud_norm{i}.Location(:, 1);
				ypt = this.ptCloud_norm{i}.Location(:, 2);

				% find closest pixel centre of each point
				x_idx = interp1(this.xc, 1:this.nx, xpt, 'nearest', 'extrap');
                y_idx = interp1(this.yc, 1:this.ny, ypt, 'nearest', 'extrap');
                this.ptCloud_norm_pixels{i} = sub2ind([this.ny, this.nx], y_idx, x_idx); % pixel indices of all points in ptCloud

			end % for cycle end

			% 3.3) and 3.4)
			for i = 1:this.lasCount
				selectedPoints = ~isnan(this.DTM(this.ptCloud_norm_pixels{i})) &... only non-NaN pixel points
								 (...
								  this.ptAttributes_norm{i}.Classification == this.VEG_LOW |... only vegetation
								  this.ptAttributes_norm{i}.Classification == this.VEG_MED |...
								  this.ptAttributes_norm{i}.Classification == this.VEG_HIGH |...
								  this.ptAttributes_norm{i}.Classification == this.GROUND... and ground points
								 );
				
				% select only vegetation and ground points from non-NaN pixels with attributes and pixel indices
				this.ptCloud_norm{i}         = select(this.ptCloud_norm{i}, selectedPoints);
				this.ptAttributes_norm{i}    = lidarPointAttributes("Classification",...
									           this.ptAttributes_norm{i}.Classification(selectedPoints));
				this.ptCloud_norm_pixels{i}  = this.ptCloud_norm_pixels{i}(selectedPoints);

			end

			time = toc(time);
			fprintf('done in %0.2f s\n', time)

		end % end of function computePointCloudAttributes

		function this = normalizePtCloud(this, options)
			arguments
				this (1,1) dataPreprocessor
				options.method (1,:) string {mustBeMember(options.method,{'DTM', 'LowestPoint'})} = 'DTM'
			end
			
			time = tic;
			
			if strcmp(options.method, "DTM")
				for i = 1:this.lasCount
% 					vegetation = ismember(this.ptAttributes{i}.Classification, [this.VEG_LOW, this.VEG_MED, this.VEG_HIGH]);
					
					X = this.ptCloud_norm{i}.Location(:, 1);
					Y = this.ptCloud_norm{i}.Location(:, 2);
					Z = this.ptCloud_norm{i}.Location(:, 3);
					
					Z = Z - this.DTM(this.ptCloud_norm_pixels{i});

					Z(Z < 0) = 0;
	
					this.ptCloud_norm{i} = pointCloud([X, Y, Z], "Intensity", this.ptCloud_norm{i}.Intensity);
				end
			end

			if strcmp(options.method, "LowestPoint")
				
				for i = 1:this.lasCount % iterate over all point clouds
% 					vegetation = ismember(this.ptAttributes{i}.Classification, [this.VEG_LOW, this.VEG_MED, this.VEG_HIGH]);
					
					X = this.ptCloud_norm{i}.Location(:, 1);
					Y = this.ptCloud_norm{i}.Location(:, 2);
					Z = this.ptCloud_norm{i}.Location(:, 3);
					% TODO: pridat "Intensity" do normalizovaneho point cloudu

					for j = 1:length(this.DTM(:)) % iterate through each pixel of mesh
						pixelPoints = this.ptCloud_norm_pixels{i}(:) == j;
						if nnz(pixelPoints) == 0
							continue;
						end
						
						Z(pixelPoints) = Z(pixelPoints) - min(Z(pixelPoints));
					end
		
					this.ptCloud_norm{i} = pointCloud([X, Y, Z], "Intensity", this.ptCloud_norm{i}.Intensity);
				end
			end
			
			time = toc(time);
			
			fprintf("Normalize points done in %0.3f s\n", time);
		end % end of normalize point cloud

		function plotPtCloud3D(this, colorMap, selectedClasses, options)
            arguments
				this (1,1) dataPreprocessor
				colorMap (:, 3)
				selectedClasses
				options.useData (1,:) string {mustBeMember(options.useData,{'original', 'normalized'})} = 'original'
			end

			if this.lasCount == 0
				error('Nothing to plot')
			end
        
			% plot only points from selected classes
			hold on

			if strcmp(options.useData, 'original')

				for i = 1:this.lasCount
					classMember = ismember(this.ptAttributes{i}.Classification, selectedClasses);
					if any(classMember)
						colorData = reshape(label2rgb(this.ptAttributes{i}.Classification, colorMap, 'k'), [], 3);
						pcshow(this.ptCloud_0{i}.Location(classMember, :), colorData(classMember, :))
					else
						warning('Nothing to plot in point cloud %d', i);
					end
				end
			
			elseif strcmp(options.useData, 'normalized')
				
				for i = 1:this.lasCount
					classMember = ismember(this.ptAttributes_norm{i}.Classification, selectedClasses);
					if any(classMember)
						colorData = reshape(label2rgb(this.ptAttributes_norm{i}.Classification, colorMap, 'k'), [], 3);
						pcshow(this.ptCloud_norm{i}.Location(classMember, :), colorData(classMember, :))
					else
						warning('Nothing to plot in point cloud %d', i);
					end
				end
			end
			
			hold off
			
		end % end of plot point cloud 3D
	end
end