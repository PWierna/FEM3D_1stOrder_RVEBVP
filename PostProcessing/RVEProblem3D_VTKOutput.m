function RVEProblem3D_VTKOutput( DisplM , StrainM , Displfluct , Strains_gips , FluctStrains_gips , Stresses_gips , IntVars_gips , MODEL , ProjectPath , Options )

%Hard settings:
dofpn = 3;

%Get mesh data:
ndofs = dofpn*size(MODEL.Coordinates,1);

%Parse fluctuating displacements:
Displfluct = Displfluct([ (1:dofpn:ndofs)' , (2:dofpn:ndofs)' , (3:dofpn:ndofs)' ]);
% if Options.insertion_type == "Displacements" %In this case, "Displfluct" has the total displacements
%     Displfluct = Displfluct - DisplM;
% end

%Perform least-squares global smoothing of variables:
nc = size(Strains_gips,2); %aux: number of strain/stress components
nodalvars = Hexa_LSGlobalSmoothing( [ Strains_gips , FluctStrains_gips , Stresses_gips , IntVars_gips ] , MODEL , Options.quadrature , 1 );
idx_nodalstrains       = 1:nc;
idx_nodalfluctstrains  = (nc+1):(2*nc);
idx_nodalstresses      = (2*nc+1):(3*nc);
idx_nodaldamage        = (3*nc+1):(3*nc+size(IntVars_gips,2));

%Assemble cell data structure:
CellData.Material_type.Values = MODEL.Conectivity(:,2);
CellData.Material_type.Type   = 'int';

%Assemble point data structure:
PointData.Displacement_Macro.Values = DisplM;
PointData.Displacement_Macro.Type   = 'float';
PointData.Displacement_Fluct.Values = Displfluct;
PointData.Displacement_Fluct.Type   = 'float';
PointData.Displacement_Total.Values = DisplM + Displfluct;
PointData.Displacement_Total.Type   = 'float';
PointData.Strain_nodal.Values = nodalvars(:,idx_nodalstrains);
PointData.Strain_nodal.Type   = 'float';
PointData.StrainFluct_nodal.Values = nodalvars(:,idx_nodalfluctstrains);
PointData.StrainFluct_nodal.Type   = 'float';
PointData.Stress_nodal.Values = nodalvars(:,idx_nodalstresses);
PointData.Stress_nodal.Type   = 'float';
PointData.Damage_nodal.Values = nodalvars(:,idx_nodaldamage);
PointData.Damage_nodal.Type   = 'float';

%Write VTK file:
Hexa_WriteVTK( MODEL.Conectivity(:,3:end) , MODEL.Coordinates(:,2:end) , CellData , PointData , ProjectPath+Options.vtk_filename );
%PointData.Points_coordinates.Values = MODEL.Coordinates(:,2:end);
%PointData.Points_coordinates.Type   = 'float';
%IPs_WriteVTK( MODEL.Coordinates(:,2:end) , PointData , Options.vtk_filename+"IPs" );
end

















