function [ Ktg , Fint , Stresses , Strains , FluctStrains , IntVarsNew ] = Hexa_StiffnessFint( U , StrainM , U_M , IntVarsPrev , RVEMODEL , MATERIALS , Options )


%% PARSE GLOBAL DATA & FIRST COMPUTATIONS ======================================== %%

%-Hard settings:
dofpn = 3; %Dofs per node

%-Get mesh data:
n_elem  = size(RVEMODEL.Conectivity,1);          %Number of Elements
npe     = size(RVEMODEL.Conectivity(:,3:end),2); %Nodes per element
dofpe   = npe*dofpn;                             %Dofs per element
ndofs   = dofpn * size(RVEMODEL.Coordinates,1);  

%-Get quadrature data:
[ n_lips , Wlips , Clips ] = QuadratureData3D( Options.quadrature );

%-Compute local interpolating data:
switch npe
    case 8
        %elem_type = "3DHexa8N";
        InterpData = LagrangeInterp3D_C08N_Local(Clips);
    case 27
        %elem_type = "3DHexa27N";
        InterpData = LagrangeInterp3D_C127N_Local(Clips);
end

%-Global IP Conectivity : (WE ASSUME ALL ELEMENTS USE THE SAME QUADRATURE)
n_gips  = n_lips * n_elem; %Total number of global ips 
GIPsCon = [repelem((1:n_elem)',n_lips,1),repmat((1:n_lips)',n_elem,1)]; %[ Element N° , Local IP N° ]

%-Compute Interpolation at all IPs:
[ ~ , detJ_gips , ~ , ~ , dNdX_gips ] = LagrangeInterp3D_Global_V( GIPsCon , InterpData , RVEMODEL );
if any(detJ_gips<=0.0)
    warning('Found non-positive Jacobian determinant');
    disp(detJ_gips(detJ_gips<=0.0));
end

%-Compute global weights for all ips:
Wgips = Wlips(GIPsCon(:,2))'.*detJ_gips;

%-Dofs matrix:
dofs_E = repelem( RVEMODEL.Conectivity(:,3:end) , 1 , dofpn )*dofpn - repmat( (dofpn:-1:1)-1 , 1 , npe ); 


%% ARRAY'S INITIALIZATIONS ======================================================= %%
%-Initialize Elemental Arrays
Fint_E = zeros( n_elem , dofpe );   %Matrix with Elemental Force Vectors 
Ktg_E  = zeros( n_elem , dofpe^2 ); %Matrix with Elemental Tangent Stiffness Matrix coeffs 


%% EVALUATE STRAINS AND CONSTITUTIVE BEHAVIOR AT ALL IPS ========================= %%

%-Compute Fluctuating Strains at (Global)IPs:
FluctStrains = ComputeCurrentStrains_at_ips( U , dNdX_gips , GIPsCon , RVEMODEL );

%-Compute Macro Strains contribution at (Global)IPs:
if Options.insertion_type == "Strains"
    macrostrains_gips = MicroBVP_InsertMacroStrains_at_ips( StrainM , GIPsCon );
elseif Options.insertion_type == "Displacements"
    UM = zeros( ndofs , 1 );
    UM(1:3:end) = U_M(:,1);    
    UM(2:3:end) = U_M(:,2);
    UM(3:3:end) = U_M(:,3);
    macrostrains_gips = ComputeCurrentStrains_at_ips( UM , dNdX_gips , GIPsCon , RVEMODEL );
end

%-Compute Total Strains at (Global)IPs:
Strains = macrostrains_gips + FluctStrains;

%-Evaluate Constitutive Response at IPs:
[ Stresses , TgModulis , IntVarsNew ] = EvaluateIndividualMaterialResponses3D( Strains , IntVarsPrev , GIPsCon , RVEMODEL , MATERIALS );


%% COMPUTE AND ASSEMBLE TG STIFFNESS MATRIX AND INTERNAL FORCE VECTOR ============ %%
for gip = 1:n_gips
    
    %-Get element:
    e = GIPsCon(gip,1);
    
    %-Get sh. functions derivatives at this ip:
    dNdX_ip   = permute(dNdX_gips(gip,:,:),[3,2,1]);
 
    %-Evaluate Strain Matrix:
    B_ip = zeros(6,dofpe);
    B_ip(1,1:dofpn:dofpe) = dNdX_ip(:,1);%dNdx
    B_ip(2,2:dofpn:dofpe) = dNdX_ip(:,2);%dNdy
    B_ip(3,3:dofpn:dofpe) = dNdX_ip(:,3);%dNdz
    B_ip(4,1:dofpn:dofpe) = dNdX_ip(:,2);%dNdy
    B_ip(4,2:dofpn:dofpe) = dNdX_ip(:,1);%dNdx
    B_ip(5,1:dofpn:dofpe) = dNdX_ip(:,3);%dNdz
    B_ip(5,3:dofpn:dofpe) = dNdX_ip(:,1);%dNdx
    B_ip(6,2:dofpn:dofpe) = dNdX_ip(:,3);%dNdz
    B_ip(6,3:dofpn:dofpe) = dNdX_ip(:,2);%dNdy
    
    %-Get tg constitutive matrix:
    Ctg_ip = reshape(TgModulis(gip,:),6,6);
    
    %-Evaluate and Weight for this (Macro) Integration Point:
    Fint_E(e,:) = Fint_E(e,:) + ( B_ip' * Stresses(gip,:)'  * Wgips(gip))'; %Internal Force
    Ktg_E(e,:)  = Ktg_E(e,:) + reshape( B_ip' * Ctg_ip * B_ip  * Wgips(gip) , 1 , dofpe*dofpe ); %Tg Stiffness

end

%Assemble Ktg indexing vectors:
idx_Ktg(:,1) = reshape( repmat(dofs_E,1,dofpe) , (dofpe^2)*n_elem , 1 );  %Row indexes
idx_Ktg(:,2) = reshape( repelem(dofs_E,1,dofpe) , (dofpe^2)*n_elem , 1 ); %Col indexes

%Assemble Tangent Stiffness Matrix and Internal Force Vector:
Fint = sparse( reshape(dofs_E,(dofpe*n_elem),1) , ones((dofpe*n_elem),1) , reshape(Fint_E,(dofpe*n_elem),1) );
Ktg  = sparse( idx_Ktg(:,1) , idx_Ktg(:,2) , reshape(Ktg_E,(dofpe*dofpe*n_elem),1) );

end



%% AUX ELEMENT FUNCTIONS ========================================================= %%
function Strains_gips = ComputeCurrentStrains_at_ips( U , dNdX_gips , GIPsCon , MODEL )

%Hard setting:
dofpn = 3;

%Parse displacements:
Ux = U(1:dofpn:end);
Uy = U(2:dofpn:end);
Uz = U(3:dofpn:end);

%Parse shape functions derivatives:
dNdx_gips = permute( dNdX_gips(:,1,:) , [1,3,2] );
dNdy_gips = permute( dNdX_gips(:,2,:) , [1,3,2] );
dNdz_gips = permute( dNdX_gips(:,3,:) , [1,3,2] );

%Get global Node numbers for each ip:
ncon_gips = MODEL.Conectivity( GIPsCon(:,1) , 3:end); 

%Compute Voigt Strain vectors for all ips:
Strains_gips = zeros( size(GIPsCon,1) , 6 );
Strains_gips(:,1:3) = [ sum(dNdx_gips.*Ux(ncon_gips),2) , sum(dNdy_gips.*Uy(ncon_gips),2) , sum(dNdz_gips.*Uz(ncon_gips),2) ]; %[xx,yy,zz]
Strains_gips(:,4)   = sum(dNdy_gips.*Ux(ncon_gips),2) + sum(dNdx_gips.*Uy(ncon_gips),2); %[xy] 
Strains_gips(:,5)   = sum(dNdz_gips.*Ux(ncon_gips),2) + sum(dNdx_gips.*Uz(ncon_gips),2); %[xz] 
Strains_gips(:,6)   = sum(dNdz_gips.*Uy(ncon_gips),2) + sum(dNdy_gips.*Uz(ncon_gips),2); %[yz] 

end

function StrainsM_gips = MicroBVP_InsertMacroStrains_at_ips( StrainM , GIPsConnect )

n_gips  = size(GIPsConnect,1);
StrainM = reshape( StrainM , 1 , length(StrainM) ); %Assert it is row-fashion

StrainsM_gips = StrainM .* ones(n_gips,1);

end


function Coordinates_gips = compute_globalipcoords_isoparam( N_gips , GIPsConnect , MODEL )

%Parse nodal coordinates to allow later indexing:
ncoords_x = MODEL.Coordinates(:,1+1);
ncoords_y = MODEL.Coordinates(:,1+2);
ncoords_z = MODEL.Coordinates(:,1+3);

%Get nodes connectivities for all ips (node number goes to 3rd dim):
ncon_gips  = permute( MODEL.Conectivity(GIPsConnect(:,1),3:end) , [1,3,2] ); 

%Compute cartesian coordinates of ips according to isoparametric interpolation:             
Coordinates_gips = sum( N_gips.*[ncoords_x(ncon_gips) ,...
                                 ncoords_y(ncon_gips) ,...
                                 ncoords_z(ncon_gips) ] , 3 ); %[ X , Y , Z ]

end


