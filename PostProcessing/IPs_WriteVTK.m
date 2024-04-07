function IPs_WriteVTK( PointsCoordinates , PointData , filename )


%% 1st COMPUTATIONS ============================================================== %%
n_points = size(PointsCoordinates,1);    %Get number of points


%% WRITE VTK FILE ================================================================ %%

%Initialize:
fid = fopen((filename + ".vtk"),'w');
fprintf(fid, '# vtk DataFile Version 4.2\n'); %In version 5.0 there is a problem w/CELL_TYPES
fprintf(fid, 'RVE MODEL_3D\n'); %HEADER
fprintf(fid, 'ASCII\n'); 
fprintf(fid, 'DATASET POLYDATA\n'); 

%Print Points coordinates:
fprintf(fid, 'POINTS %d float\n', n_points); 
fprintf(fid, [repmat(number_format('float'),1,3),' \r\n'], PointsCoordinates'); %Remember to transpose when printing

%Point Data
if isstruct( PointData ) %If PointData is not empty[]
    
    %Get Point Data:
    point_data_fnames  = fieldnames( PointData );   %Get field names
    point_data_nfields = length(point_data_fnames);

    %Print fields for point data:
    fprintf(fid, '\nPOINT_DATA %d \n', n_points);
    fprintf(fid, 'FIELD FieldData %d \n',point_data_nfields); 
    for f = 1:point_data_nfields %(For each field in PointData struct)
        field_name   = point_data_fnames{f};
        field_values = PointData.(field_name).Values;
        field_type   = PointData.(field_name).Type;
        field_comps  = size(field_values,2); %field components
        fprintf(fid, [ field_name , ' %d %d ' , field_type , '\n' ] , field_comps , n_points); 
        fprintf(fid, [ repmat(number_format(field_type),1,field_comps),' \r\n'], field_values');
    end
    
end

%Close file
fclose(fid);

end


function nformat = number_format( field_type )

switch field_type
    case 'float'
        nformat = '%10.6g ';
    case 'int'
        nformat = '%d ';
end

end