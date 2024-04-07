function MacroDisplacements = MicroBVP3D_InsertMacroDisplacements_at_nodes( DisplM , StrainM , RVEMODEL )

%Modify Macro-Strain vector from engineering to tensorial:
StrainM(4:6) = StrainM(4:6)./2;  

%Get nodal coordinates:
ncoords  = RVEMODEL.Coordinates(:,2:end);

%Assumming prismatic rve, compute Delta_ncoords (relative to the geometric centre):
Dncoords = ncoords - 1/2*(max(ncoords)+min(ncoords)); %Only valid for prismatic RVEs!!

%Compute Components of the Macro Displacements:
UMx = DisplM(1) + StrainM(1) * Dncoords(:,1) + StrainM(4) * Dncoords(:,2)  + StrainM(5) * Dncoords(:,3) ;
UMy = DisplM(2) + StrainM(4) * Dncoords(:,1) + StrainM(2) * Dncoords(:,2)  + StrainM(6) * Dncoords(:,3) ;
UMz = DisplM(3) + StrainM(5) * Dncoords(:,1) + StrainM(6) * Dncoords(:,2)  + StrainM(3) * Dncoords(:,3) ;
  
%Assemble output:  
MacroDisplacements = [ UMx , UMy , UMz ];

end