function parameter_values = coord_op_param_from_EPSG_database(EPSG_coordinate_operation_code)

%% Make connection to database
try 
    conn = database('EPSG database','','');
% 	conn = database('data/EPSG-v10_033-Access','','');
catch
    configureODBCDataSource()
    error('Add ''EPSG database'' using ODBC Data Source Administrator. Use Data Source name ''EPSG database'' and select file ''EPSG-v10_033-Access.mdb'' in your data directory.')
end

%% Execute query and fetch results
% Coordinate_Operation = fetch(conn,'SELECT * FROM `Coordinate_Operation`');
% Coordinate_Operation_Method = fetch(conn,'SELECT * FROM `Coordinate_Operation Method`');
% Coordinate_Operation_Parameter = fetch(conn,'SELECT * FROM `Coordinate_Operation Parameter`');
% Coordinate_Operation_Parameter_Usage = fetch(conn,'SELECT * FROM `Coordinate_Operation Parameter Usage`');
Coordinate_Operation_Parameter_Value = fetch(conn,'SELECT * FROM `Coordinate_Operation Parameter Value`');

%% Close connection to database
close(conn)

%% Get coordinate operation parameter values
ind_parameter_value = find(Coordinate_Operation_Parameter_Value.COORD_OP_CODE == EPSG_coordinate_operation_code);
parameter_values = Coordinate_Operation_Parameter_Value.PARAMETER_VALUE(ind_parameter_value);

% coordinate_operation_code = 8365;

% ind_op = find(Coordinate_Operation.COORD_OP_CODE == coordinate_operation_code);
% Coordinate_Operation(ind_op,:)

% method_code = Coordinate_Operation.COORD_OP_METHOD_CODE(ind_op);
% ind_method = find(Coordinate_Operation_Method.COORD_OP_METHOD_CODE == method_code);
% Coordinate_Operation_Method(ind_method,:)

% Coordinate_Operation_Parameter_Value(ind_parameter_value,:)
% parameter_code = Coordinate_Operation_Parameter_Value.PARAMETER_CODE(ind_parameter_value);

end
