function varargout = JTSK03_to_gps_transformation(x, y, h)
	if (nargin == 2)
		h = 0;
	end
	
	% S-JTSK(JTSK03) - Krovak East North (x,y) to Geodetic S-JTSK(JTSK03) (lat, lon, h)
	JTSK03 = projcrs(8353);
	[lat, lon] = projinv(JTSK03, x, y);

	% Geodetic S-JTSK(JTSK03) (lat, lon, h) to Geocentric S-JTSK(JTSK03) (X,Y,Z)
	[X, Y, Z] = geodetic2ecef(referenceEllipsoid('Bessel 1841'), lat, lon, h);
	
	% datum transformation: S-JTSK(JTSK03) to WGS84
% 	code = 5239; % EPSG code of coordinate operation
%     p = coord_op_param_from_EPSG_database(code);
%     p(4:6) = p(4:6)*(pi/180/60/60); % conversion from arc seconds to radians
%     p(7) = 1+p(7)*1e-6; % scale factor
%     XYZ = d3trafo([X',Y',Z'],p(1:7));
%     X = XYZ(:,1)';
%     Y = XYZ(:,2)';
%     Z = XYZ(:,3)';

	%#################################################################################%
	% skusal som porovnat rozdiely medzi prvym a poslednym bodom krivky,
	% kedze tie by sa nemali hybat pocas evolucie
	% pri pouziti jednej EPSG:5239 boli rozdiely medzi:
	%	- prvymi bodmi:		1.0954e-05
	%	- poslednymi bodmi: 1.0954e-05
	%
	% pri pouziti 2 transformacii boli rozdiely medzi:
	%	- prvymi bodmi:		1.0517e-07
	%	- poslednymi bodmi: 1.0514e-07
	% 
	% cize asi bude lepsie pouzit tie 2 ako v povodnej funkcii na
	% transformaciu suradnic, aj ked teda tie rozdiely su v oboch pripadoch
	% dost male
	
	% datum transformation: S-JTSK(JTSK03) to ETRS89
	code = 8365; % EPSG code of coordinate operation
    p = coord_op_param_from_EPSG_database(code);
	p = -p;
    p(4:6) = p(4:6)*(pi/180/60/60); % conversion from arc seconds to radians
    p(7) = 1+p(7)*1e-6; % scale factor
    XYZ = d3trafo([X',Y',Z'],p(1:7));
    X = XYZ(:,1)';
    Y = XYZ(:,2)';
    Z = XYZ(:,3)';
	
	% datum transformation: ETRS89 to WGS84
	code = 9225; % EPSG code of coordinate operation
    p = coord_op_param_from_EPSG_database(code);
	p = -p;
    p(4:6) = p(4:6)*(pi/180/60/60); % conversion from arc seconds to radians
    p(7) = 1+p(7)*1e-6; % scale factor
    XYZ = d3trafo([X',Y',Z'],p(1:7));
    X = XYZ(:,1)';
    Y = XYZ(:,2)';
    Z = XYZ(:,3)';
	%#################################################################################%
	
	
	% (X,Y,Z) to (lat, lon, h) in WGS84
	[lat, lon, h] = ecef2geodetic(referenceEllipsoid('wgs84'), X, Y, Z);
	
	varargout{1} = lat;
	varargout{2} = lon;
	varargout{3} = h;
end