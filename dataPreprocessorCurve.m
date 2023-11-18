classdef dataPreprocessorCurve
	% Class for preprocessing point cloud data:
	% 1) load point cloud data ptCloud_0 and curve as polyshape
	% 2) filter point cloud ptCloud_0:
	%   2.1) select points within Omega bounds
	%   2.2) select only vegetation and ground points
	% 3) normalize point cloud PC_0 -> PC_n

	properties
		% 2D mesh parameters
		x1, x2, y1, y2                  % x/y pixel centers range
        nx, ny                          % number of pixels in x/y direction
        h	                            % pixel sizes
        x, y                            % pixel edge coordinates
        xc, yc, Xc, Yc                  % pixel center coordinates (and meshgrids)
		Omega polyshape					% polyshape object representing mesh

		curve polyshape
        
		% PC stuff
		ptCloud_0						% cell of all point clouds
		ptCloud_norm				    % cell of normalized point cloud
		lasCount						% number of point clouds
		
		ptAttributes_0					% cell of classification for all points clouds
		ptAttributes_norm

		ptCloud_0_pixels				% DTM pixel indices of all points in ptCloud
		ptCloud_norm_pixels

		GROUND = 2
		VEG_LOW = 3
		VEG_MED = 4
		VEG_HIGH = 5

	end % end of properties

	methods
		% CONSTRUCTOR
		function this = dataPreprocessorCurve(curve, h, n)
			arguments
				curve polyshape
				h
				n
			end

			this.curve = curve;

			% MESH
			[xlim_c, ylim_c] = curve.boundingbox;
			BB_width  = abs(xlim_c(2) - xlim_c(1)); % sirka bounding boxu kml krivky
			BB_height = abs(ylim_c(2) - ylim_c(1)); % vyska bounding boxu kml krivky

			% vypocet sirky a vysky pre \Omega
			Omega_width  = ceil(BB_width / h) * h + n * h;
			Omega_height = ceil(BB_height / h) * h + n * h;

			% offset v smeroch x a y -> sluzi na vypocet suradnic pre vrcholy \Omega
			offset_x = (Omega_width - BB_width) / 2;
			offset_y = (Omega_height - BB_height) / 2;

			% vypocet suradnic pre vrcholy \Omega
			x_new = [xlim_c(1) - offset_x, xlim_c(2) + offset_x, xlim_c(2) + offset_x, xlim_c(1) - offset_x];
			y_new = [ylim_c(2) + offset_y, ylim_c(2) + offset_y, ylim_c(1) - offset_y, ylim_c(1) - offset_y];

			this.Omega = polyshape(x_new, y_new);
			[xlim_o, ylim_o] = this.Omega.boundingbox;

			% vytvorenie vypoctovej oblasti
			this.x1 = xlim_o(1);
			this.x2 = xlim_o(2);
			this.y1 = ylim_o(1);
			this.y2 = ylim_o(2);

			this.nx = Omega_width;%/ 1;
			this.ny = Omega_height;% / 1;

			this.h = 1;

			this.xc = linspace(this.x1 + this.h/2, this.x2 - this.h/2, this.nx);
			this.yc = linspace(this.y1 + this.h/2, this.y2 - this.h/2, this.ny);
			[this.Xc, this.Yc] = meshgrid(this.xc, this.yc);

			this.x = this.x1:this.h:this.x2;
			this.y = this.y1:this.h:this.y2;


		end % end of constructor
		
		function this = filterPointCloud(this, ptCloud, ptAttributes)
			% kontrola, ci je priradeny nejaky point cloud
			if isempty(this.ptCloud_0)
				error('No point cloud available');
			end

			fprintf('Filtering point cloud ... ');
            time = tic;

			this.ptCloud_0 = ptCloud;
			this.ptAttributes_0 = ptAttributes;

			% 2.1), 2.2)
			selectedPoints = this.ptCloud_0.Location(:,1) >= this.x1 &... only points within Omega
							 this.ptCloud_0.Location(:,1) <= this.x2 &...
							 this.ptCloud_0.Location(:,2) >= this.y1 &...
							 this.ptCloud_0.Location(:,2) <= this.y2 &...
							 ( ... only vegetation and ground points
							 this.ptAttributes_0.Classification == this.VEG_LOW  |... only vegetation
							 this.ptAttributes_0.Classification == this.VEG_MED  |...
							 this.ptAttributes_0.Classification == this.VEG_HIGH |...
							 this.ptAttributes_0.Classification == this.GROUND ... and ground points
							 );

			this.ptCloud_norm      = select(this.ptCloud_0, selectedPoints);
			this.ptAttributes_norm = lidarPointAttributes("Classification",...
														  this.ptAttributes_0.Classification(selectedPoints));
			
			time = toc(time);
			fprintf('done in %0.2f s\n', time)

		end % end of function computePointCloudAttributes

		function this = normalizePtCloud(this)
			arguments
				this (1,1) dataPreprocessorCurve
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

		function plotMesh(this, options)
			arguments
				this dataPreprocessorCurve
				options
			end

			hold on
			plot(this.curve, "EdgeColor", "blue", "FaceAlpha", 0.05, "FaceColor","blue")
			plot(this.Omega,"FaceAlpha",0,"LineStyle","-","LineWidth",2)
			for i = 1:length(this.x)
				plot([this.x(i) this.x(i)],[this.y(1) this.y(end)],...
					'Color', "#EDB120", 'LineWidth', 0.5,'LineStyle','-') %y grid lines
			end

			for i = 1:length(this.y)
				plot([this.x(1) this.x(end)],[this.y(i) this.y(i)],...
					'Color', "#EDB120", 'LineWidth', 0.5,'LineStyle','-') %x grid lines
			end

			x = this.curve.Vertices(:,1);
			y = this.curve.Vertices(:,2);
			[XYin, ~] = inpolygon(this.Xc, this.Yc, x , y);

			scatter(this.Xc(XYin),  this.Yc(XYin), 5, "*g")
			scatter(this.Xc(~XYin), this.Yc(~XYin),5, "*r")
			
			hold off
			axis equal
		end

		function plotPtCloud3D(this, colorMap, selectedClasses, options)
            arguments
				this (1,1) dataPreprocessorCurve
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
					classMember = ismember(this.ptAttributes_0{i}.Classification, selectedClasses);
					if any(classMember)
						colorData = reshape(label2rgb(this.ptAttributes_0{i}.Classification, colorMap, 'k'), [], 3);
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