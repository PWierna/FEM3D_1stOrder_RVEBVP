function Hexa_WriteVTK( ModelConnectivities , PointsCoordinates , CellData , PointData , filename )


%% 1st COMPUTATIONS ============================================================== %%
n_cells  = size(ModelConnectivities,1);  %Get number of cells
n_points = size(PointsCoordinates,1);    %Get number of points
npe = size(ModelConnectivities,2);       %Points per cell
dofpn = 3;                               %Number of dofs per node


%% ASSEMBLE POINT AND CELL MATRICES ============================================== %%

%Resort connectivities and assign cell types:
switch npe
    case 8  %Hexa8N
        sort_connect = 1:8;
        cell_type = 12 * ones(1,n_cells);         %12 = n° cell type for 8Nhexahedron in VTK
    case 27 %Hexa27N
        sort_connect = [1:8,9:12,17:20,13:16,25,23,22,24,21,26,27];
        cell_type = 29 * ones(1,n_cells);         %29 = n° cell type for 27Nhexahedron in VTK
end

%Assemble Cells Connectivities:
CellsConnectivities = [ npe*ones(n_cells,1) (ModelConnectivities(:,sort_connect)-1) ]; %(1st col=N°points defining each cell)

%Get total N° of values defining the cells:
cell_numb = numel(CellsConnectivities); 


%% WRITE VTK FILE ================================================================ %%

%Initialize:
fid = fopen((filename + ".vtk"),'w');
fprintf(fid, '# vtk DataFile Version 4.2\n'); %In version 5.0 there is a problem w/CELL_TYPES
fprintf(fid, 'RVE MODEL_3D\n'); %HEADER
fprintf(fid, 'ASCII\n'); 
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n'); 

%Print Points coordinates:
fprintf(fid, 'POINTS %d float\n', n_points); 
fprintf(fid, [repmat(number_format('float'),1,dofpn),' \r\n'], PointsCoordinates'); %Remember to transpose when printing

%Print Cells connectivities:
fprintf(fid, '\nCELLS %d %d\n', n_cells , cell_numb); 
fprintf(fid, [repmat('%d ',1,npe+1),' \r\n'], CellsConnectivities'); 
fprintf(fid, '\nCELL_TYPES %d \n', n_cells); 
fprintf(fid, '%d \r\n', cell_type);

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
 
%Cell Data
if isstruct( CellData )
        
    %Get Cell Data:
    cell_data_fnames  = fieldnames( CellData );   %Field names
    cell_data_nfields = length(cell_data_fnames); %Number of fields
    
    %Print fields for cell data:
    fprintf(fid, '\nCELL_DATA %d \n', n_cells);
    fprintf(fid, 'FIELD FieldData %d \n',cell_data_nfields);
    for f = 1:cell_data_nfields %(For each field in CellData struct)
        field_name   = cell_data_fnames{f};
        field_values = CellData.(field_name).Values;
        field_type   = CellData.(field_name).Type;
        field_comps  = size(field_values,2); %field components
        fprintf(fid, [ field_name , ' %d %d ' , field_type , '\n' ] , field_comps , n_cells); 
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