function [ NodalValues ] = Hexa_LSGlobalSmoothing( IPValues , MODEL , Quadrature , UseLumped )


%% PARSE GLOBAL DATA & FIRST COMPUTATIONS ======================================== %%

%-Parse input data:
n_nodes = size(MODEL.Coordinates,1);             %Number of nodes
n_elem  = size(MODEL.Conectivity,1);             %Number of elements
npe     = size(MODEL.Conectivity(:,3:end),2);    %Nodes per element
ncomps  = size(IPValues,2);                      %Components to be smoothed

%-Get quadrature data:
[ n_lips , Wlips , Clips ] = QuadratureData3D( Quadrature );

%-Compute local interpolating data:
switch npe
    case 8 %"Hexa8N";
        InterpData = LagrangeInterp3D_C08N_Local(Clips);
    case 27%"Hexa27N";
        InterpData = LagrangeInterp3D_C127N_Local(Clips);
end

%-Global IP Conectivity : (WE ASSUME ALL ELEMENTS USE THE SAME QUADRATURE)
%n_gips  = n_lips * n_elem; %Total number of global ips 
GIPsCon = [repelem((1:n_elem)',n_lips,1),repmat((1:n_lips)',n_elem,1)]; %[ Element N° , Local IP N° ]

%-Compute Interpolation at all IPs:
[ N_gips , detJ_gips , ~ , ~ , ~ ] = LagrangeInterp3D_Global_V( GIPsCon , InterpData , MODEL );

%-Compute global weights for all ips:
Wgips = Wlips(GIPsCon(:,2))'.*detJ_gips;

%-Initialize Global Entities:
S       = sparse( n_nodes , n_nodes );              %Smoothing Matrix
S_E     = zeros( npe , npe , n_elem );              %Elemental Smoothing Matrices
SVect   = zeros( n_nodes , size(IPValues,2) );      %Smoothing vector's matrix 
SVect_E = zeros( npe , ncomps , n_elem );           %Elemental Smoothing Vectors
gip   = 0;

%% LOOP OVER ELEMENTS ============================================================ %%      
for e = 1:n_elem 
    
    %-Get Element Nodes, Coordinates & DoFs:
    con_e = MODEL.Conectivity(e,3:end);     %Global Node numbers
    
    %-Initialize elemental entities
    S_e = zeros( npe );
    SVect_e = zeros( npe , size(IPValues,2) );
    
    %--[Loop over Integration Points] --->----------------------------------------->%
    for ip = 1:n_lips

        %-Global IP Counter:
        gip = gip + 1;
        
        %-Get sh. functions and the Jacobian determinant at this ip:
        N_ip = permute(N_gips(gip,:,:),[1,3,2]);
        
        %-Integrate Elemental Smoothing Matrix and Smoothing "Force" Vectors:   
        S_e     = S_e + N_ip' * N_ip * Wgips(gip);
        SVect_e = SVect_e + N_ip' * IPValues(gip,:) * Wgips(gip);
               
        %Accumulate:
        S_E(:,:,e) = S_E(:,:,e) + N_ip' * N_ip  * Wgips(gip);
        SVect_E(:,:,e) = SVect_E(:,:,e) + N_ip' * IPValues(gip,:) * Wgips(gip);
        
    end %<--------------------------------<---[Ends Loop over Integration Points] --%
    
    %-Assemble Smoothing Matrix and Smoothing Vectors:
    S( con_e , con_e ) = S( con_e , con_e ) + S_e;
    SVect( con_e , : ) = SVect( con_e , : ) + SVect_e;
    
end %<--------------------------------------------- <--- [Ends Loop over Elements]--%

%Assembly:
Sidx_cols = repmat( permute(MODEL.Conectivity(:,3:end),[3,2,1]) , npe , 1 , 1 );
Sidx_rows = repmat( permute(MODEL.Conectivity(:,3:end),[2,3,1]) , 1 , npe , 1 );

SVectidx_rows = repmat( permute(MODEL.Conectivity(:,3:end),[2,3,1]) , 1 , ncomps , 1 );
SVectidx_cols = repmat( (1:ncomps) , npe , 1  , n_elem );

S2 = sparse( reshape(Sidx_cols,npe*npe*n_elem,1) ,...
             reshape(Sidx_rows,npe*npe*n_elem,1) ,...
             reshape(S_E,npe*npe*n_elem,1) );
         
SVect2 = sparse( reshape(SVectidx_cols,npe*ncomps*n_elem,1) ,...
                 reshape(SVectidx_rows,npe*ncomps*n_elem,1) ,...
                 reshape(SVect_E,npe*ncomps*n_elem,1) );

             
%% COMPUTE LEAST-SQUARES GLOBAL SMOOTHING ======================================== %%
if UseLumped == 1
    SLump = sum(S,2); %Lumped Smoothing matrix
    NodalValues  = (1./SLump) .* SVect ;
else
    NodalValues = S \ SVect ;
end

return