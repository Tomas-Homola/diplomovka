classdef curve
	% class of dictrete curves
	properties
		open                     % 0 = closed curve, 1 = open curve
		xy                       % coordinates of points
		Npts                     % number of points
		iprev,inext,iprev2,inext2  % array for previous and next indices
		r                        % vector to the next point
		r2                       % vector from previous to next point
		h                        % edge length
		L                        % length of a curve
		area                     % area (of simple closed curve)
		T,N                      % tangent and normal vectors
		ksig                     % signed curvature at curve points
		ksigEdge                 % signed curvature at midpoints of edges
		curvePixels              % pixels occupied by curve points
	end

	methods
		function C = curve(topology)
			% constructor
			if nargin == 0
				C.open = 0;      % default topology = closed curve
			else
				C.open = topology;
			end
		end

		function C = set_coordinates(C,xy)
			% set (x,y) coordinates of points of discrete curve C
			C.xy = xy;
			C.Npts = size(C.xy,1);
			C.iprev = circshift(1:C.Npts,1)';
			C.inext = circshift(1:C.Npts,-1)';
			C.iprev2=circshift(1:C.Npts,2);
			C.inext2=circshift(1:C.Npts,-2);
		end

		function C = circle(C,S,r,h)
			% create circle with radius r centered at S
			C.open = 0;
			n = floor((2*pi*r)/(0.5*h));
			uu = linspace(0,-2*pi,n+1)';
			u = uu(1:n);
			xy_circ = [S(1) + r*cos(u),S(2) + r*sin(u)];
			C = set_coordinates(C,xy_circ);
			C = compute_geometry(C);
		end

		function C = compute_geometry(C)
			% compute all geometric properties of curve C
			C.Npts = size(C.xy,1);
			C.r = C.xy(C.inext,:) - C.xy;                % vector to the next point
			if C.open
				C.r(C.Npts,:) = [0,0];
			end
			C.h = vecnorm(C.r,2,2);                      % edge length
			C.L = sum(C.h);                              % curve length
			C.r2 = C.xy(C.inext,:) - C.xy(C.iprev,:);    % vector from previous to next point
			if C.open
				C.r2(1,:) = 2*C.r(1,:);
				C.r2(C.Npts,:) = 2*C.r(C.Npts-1,:);
			end
			C.T = C.r2./vecnorm(C.r2,2,2);               % tangent vector
			C.N = [-C.T(:,2),C.T(:,1)];                  % positively oriented normal

			% signed curvature at midpoints of edges
			C.ksigEdge = dot(C.xy(C.iprev,:) - C.xy - C.xy(C.inext,:) + C.xy(C.inext2,:),[-C.r(:,2),C.r(:,1)],2)./(2*C.h.^3);
			if C.open
				C.ksigEdge(1) = 0;
				C.ksigEdge(C.Npts-1) = 0;
				C.ksigEdge(C.Npts) = 0;
			end

			% signed curvature at curve points
			%C.ksig = dot(4*(C.xy(C.iprev,:) - 2*C.xy + C.xy(C.inext,:))./dot(C.r2,C.r2,2),C.N,2);
			C.ksig = (C.ksigEdge(C.iprev,:)+C.ksigEdge)/2;
		end

		function C = compute_some_geometry(C)
			% recompute necessary geometric properties of curve C
			C.r = C.xy(C.inext,:) - C.xy;                % vector to the next point
			if C.open
				C.r(C.Npts,:) = [0,0];
			end
			C.h = vecnorm(C.r,2,2);                      % edge length
		end

		function C = compute_curvePixels(C,terrain)
			% find pixels of in plane occupied by the curve points
			ii = interp1(terrain.yc,1:terrain.ny,C.xy(:,2),'nearest','extrap');
			jj = interp1(terrain.xc,1:terrain.nx,C.xy(:,1),'nearest','extrap');
			C.curvePixels = sub2ind([terrain.ny,terrain.nx],ii,jj);
		end

		function C =calcArea(C)
			% TOMEK: určite to dobre funguje pre
			% jednoduchú uzavretú krivku). Bolo by fajn zistiť, ako polyarea funguje
			% pre nie jednoduché uzavreté (také zauzlené) krivky, ktoré nám
			% tiež pri splittingu vznikali.
			C.area=polyarea(C.xy(:,1),C.xy(:,2));
		end

		function plot(C, options)
			arguments
				C curve
				options.Color = 'red'
				options.MarkerSize = 10;
			end
			% plot curve

			if C.open
				plot(C.xy(:,1),C.xy(:,2),'.-',...
					'Color',options.Color,'MarkerSize', options.MarkerSize); % povodne 10
			else
				plot([C.xy(:,1); C.xy(1,1)],[C.xy(:,2); C.xy(1,2)],'.-',...
					'Color',options.Color,'MarkerSize', options.MarkerSize);
			end
			axis equal
			%row=dataTipTextRow('index',C.xy(:,1))
			%            row = dataTipTextRow('index',[1:20]);
			%            p.DataTipTemplate.DataTipRows(end+1) =row;
		end

		function plot_T(C,scale,color)
			% plot tangent vectors
			if nargin == 1
				scale = 1;
			end
			quiver(C.xy(:,1),C.xy(:,2),scale.*C.T(:,1),scale.*C.T(:,2),color,'AutoScale','off')
		end

		function plot_N(C,scale,color)
			% plot normal vectors
			if nargin == 1
				scale = 1;
			end
			quiver(C.xy(:,1),C.xy(:,2),scale.*C.N(:,1),scale.*C.N(:,2),color,'AutoScale','off')
		end
	end
end