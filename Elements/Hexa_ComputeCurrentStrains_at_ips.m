function Strains_gips = Hexa_ComputeCurrentStrains_at_ips( U , dNdX_gips , GIPsCon , MODEL )

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

return