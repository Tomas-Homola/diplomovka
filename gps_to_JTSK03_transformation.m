function varargout = gps_to_JTSK03_transformation(lat,lon,h)
% Transformation from GPS coordinates (WGS84) to S-JTSK [JTSK03]

    if nargin == 2
        h = 0;
    end
    % GPS to XYZ on WGS 84 Ellipsoid
    [X,Y,Z] = geodetic2ecef(referenceEllipsoid('wgs84'),lat,lon,h);

    % datum transformation 1: WGS 84 to ETRS89
    code = 9225; % EPSG code of coordinate operation
    p = coord_op_param_from_EPSG_database(code);
    p(4:6) = p(4:6)*(pi/180/60/60); % conversion from arc seconds to radians
    p(7) = 1+p(7)*1e-6; % scale factor
    XYZ = d3trafo([X',Y',Z'],p(1:7));
    X = XYZ(:,1)';
    Y = XYZ(:,2)';
    Z = XYZ(:,3)';

    % datum transformation 2: ETRS89 to S-JTSK [JTSK03]
    code = 8365; % EPSG code of coordinate operation
    p = coord_op_param_from_EPSG_database(code);
    p(4:6) = p(4:6)*(pi/180/60/60); % conversion from arc seconds to radians
    p(7) = 1+p(7)*1e-6; % scale factor
    XYZ = d3trafo([X',Y',Z'],p);
    X = XYZ(:,1)';
    Y = XYZ(:,2)';
    Z = XYZ(:,3)';

    % S-JTSK [JTSK03] XYZ to latlon (geographic 2D)
    [lat1,lon1,h] = ecef2geodetic(referenceEllipsoid('Bessel 1841'),X,Y,Z);

    % S-JTSK [JTSK03] (geographic 2D) to S-JTSK [JTSK03] / Krovak East North
    JTSK03 = projcrs(8353);
    [x,y] = projfwd(JTSK03,lat1,lon1);
    
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = h;
end