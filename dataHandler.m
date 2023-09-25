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
% 			this.alphaData_DTM = ones(size(this.DTM)); % vytvorenie pola, kde budu ulozene udaje o alpha pre kazdy pixel
% 			this.alphaData_DTM(isnan(this.DTM)) = 0; % tam, kde su NaN hodnoty, tak budu priehladne

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

			time = toc(time);
			fprintf('done in %0.2f s\n', time)

		end % end of function computePointCloudAttributes

		function this = normalizePtCloud(this)
			
			this.ptCloud_normalized = cell(this.lasCount, 1);

			for i = 1:this.lasCount
				vegetation = ismember(this.ptAttributes{i}.Classification, [this.VEG_LOW, this.VEG_MED, this.VEG_HIGH]);
				
				X = this.ptCloud{i}.Location(:, 1);
				Y = this.ptCloud{i}.Location(:, 2);
				Z = this.ptCloud{i}.Location(:, 3);
				
				Z(vegetation) = Z(vegetation) - this.DTM(this.ptCloud_pixels{i}(vegetation));

				this.ptCloud_normalized{i} = pointCloud([X, Y, Z]);
			end

			fprintf("Normalize points done\n");
		end % end of normalize point cloud

		function plotPtCloud2D(this, colorData, selectedClass) 
% 			arguments
% 				this dataHandler
% 				colorData (:, 3) {mustBeNumeric}
% 				selectedClass 
% 			end
% 			
% 			
% 			if this.lasCount == 0
% 				warning('No Point Cloud')
% 				return
% 			end
% 			
% % 			figure
% 			% TH: pri vykreslovani vsetkych bodov nie je funkcia scatter zrovna najlepsia, co sa tyka plynulosti 
% % 			for i = 1:lasCount
% % 				scatter(this.ptCloud{i}.Location(:, 1), this.ptCloud{i}.Location(:, 2), 1, colorData{i}/255)
% % 				hold on
% % 			end
% 			% TH: preto sa vykresluje zatial vzdy iba jedna vybrana trieda
% 			for i = 1:this.lasCount
% 				classMember = this.ptAttributes{i}.Classification == selectedClass;
% 				if any(classMember)
% 					scatter(this.ptCloud{i}.Location(classMember, 1), this.ptCloud{i}.Location(classMember, 2),...
% 						10, colorData{i}(classMember, :)/255, 'filled')
% 					%hold on
% 				else
% 					%warning('Nothing to plot')
% 				end
% 			end
% 			
% 			%hold off
%             %title 'Terrain Point Cloud'
% %             axis ij
%             %axis equal
		end % end of plot point cloud 2D

		function plotPtCloud3D(this, colorData, selectedClasses, options)
            arguments
				this (1,1) dataHandler
				colorData
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
						pcshow(this.ptCloud{i}.Location(classMember, :), colorData{i}(classMember, :))
					else
						warning('Nothing to plot')
					end
				end
			
			elseif strcmp(options.useData, 'normalized')
				
					for i = 1:this.lasCount
					classMember = ismember(this.ptAttributes{i}.Classification, selectedClasses);
					if any(classMember)
						pcshow(this.ptCloud_normalized{i}.Location(classMember, :), colorData{i}(classMember, :))
					else
						warning('Nothing to plot')
					end
				end


			end
			
			hold off
			
		end % end of plot point cloud 3D

	end % end of methods
end







