function H = Hexa_MKCmat( MODEL , faces_connect , faces_normals , Options )


%% PARSE INPUT DATA & FIRST COMPUTATIONS ========================================= %%

%-Hard settings:
dofpn   = 3; %Dofs per node

%-Get mesh data:
npe       = size(MODEL.Conectivity(:,3:end),2); %Nodes per element
dofpe     = npe*dofpn;                          %Dofs per element
n_elem    = size(MODEL.Conectivity,1);          %Number of elements

switch npe
    case 27
        elem_type = "3DHexa27N";%element type
        %npf       = 9;          %nodes per face
    case 8
        elem_type = "3DHexa8N"; %element type
        %npf       = 4;          %nodes per face
end

%-Get equivalent 2D quadrature data:
switch Options.quadrature
    case "Gauss 3x3x3"
        [ n_ips2D , Wips2D , Cips2D ] = QuadratureData2D( "Gauss 3x3" );
    case "Gauss 2x2x2"
        [ n_ips2D , Wips2D , Cips2D ] = QuadratureData2D( "Gauss 2x2" );
end

%-Compute dummy local interpolating data 3D to get faces local connectivities:
switch elem_type
    case "3DHexa8N"
        %TODO
    case "3DHexa27N"
        [ ~ , ~ , dummyCips3D ] = QuadratureData3D( Options.quadrature );
        [ ~ , LFacesConnect ]   = LagrangeInterp3D_C127N_Local(dummyCips3D);
end

%-Get Face Elements connectivities and normals:
% faces_connect = [ 1*ones(size(MODEL.BoundingFaces.XNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.XNormal.Minus.Connectivities ;
%                   2*ones(size(MODEL.BoundingFaces.XNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.XNormal.Plus.Connectivities ];
% faces_normals = [ MODEL.BoundingFaces.XNormal.Minus.Normal' ;
%                   MODEL.BoundingFaces.XNormal.Plus.Normal'  ];
n_facelem     = size(faces_connect,1);

%-Initialize DOF Arrays:
%H = zeros( 6 , dofpn*size(MODEL.Coordinates,1) );   %H Matrix ( F I X M E !!! )
H_E    = zeros( 6 , dofpe , n_elem );

%-Dofs matrix:
dofs_E = repelem( MODEL.Conectivity(:,3:end) , 1 , dofpn )*dofpn - repmat( [dofpn:-1:1]-1 , 1 , npe ); 

%-Resort index for local nodal pos at faces:
faceidxs_resort2D = Hexa27Nresortfacenodes2Quad9N;
    

%% LOOP OVER ELEMENTAL FACES AT THE BOUNDARIES =================================== %%
for ef = 1:n_facelem
    
    %Get the element this face belongs to:
    e = faces_connect(ef,2);
    
    %Identify corresponding local face:
    lface_idx = find( all( LFacesConnect.Connectivities == faces_connect(ef,3:end) , 2 ) );
    
    %Get (global/cartesian) normal to this face:
    ef_normal = faces_normals( faces_connect(ef,1) , : );
    
    %Assemble n matrix for this face (assumed constant over the whole face):
    ef_nmat   = [             diag(ef_normal)            ;
                  ef_normal(2) ef_normal(1)      0       ;
                  ef_normal(3)       0      ef_normal(1) ;
                        0      ef_normal(3) ef_normal(2) ];
    
    %Assemble 3D Coordinates for 2D Gauss Quadrature over the face:
    lface_ipcoords = zeros(n_ips2D,3);
    lface_ipcoords(:,LFacesConnect.Normals(lface_idx,:)~=0) = LFacesConnect.Normals( lface_idx , LFacesConnect.Normals(lface_idx,:)~=0 );
    lface_ipcoords(:,LFacesConnect.Normals(lface_idx,:)==0) = Cips2D;
    
    %Get shape functions values for all the IPs at this face:
    switch elem_type
        case "3DHexa8N"
            %TODO
        case "3DHexa27N"
            [ interpdata3D , ~ ] = LagrangeInterp3D_C127N_Local(lface_ipcoords);
            N3D_ips = permute( interpdata3D(:,1,:) , [3,1,2] );
    end

    %Initialize elemental H matrix:
    He = zeros(6,dofpe);
    
    
    %Perform integration over the face -------------------------------------------- %
    
    %-Get equivalent 2D element for this face:
    ef_2Dnodes  = MODEL.Conectivity( e , 2 + faces_connect(ef,3:end) );
    ef_2Dnodes  = ef_2Dnodes( faceidxs_resort2D(lface_idx,:) ); %Resort
    ef_2Dcoords = MODEL.Coordinates( ef_2Dnodes , 1+find(ef_normal==0) );
    
    %-Run check-loop over ips:
    flag_flip = 0;
    for ip2D = 1:n_ips2D
        [ ~ , ~ , detJ_ip2D , ~ , ~ ] = LagrangeInterp2D_C19N( ef_2Dcoords , Cips2D(ip2D,:) );
        if detJ_ip2D <= 0; flag_flip=1; break
        end
    end
    
    if flag_flip == 1;
        ef_2Dnodes = ef_2Dnodes([1 8:-1:2 9]); %flip order
        ef_2Dcoords = MODEL.Coordinates( ef_2Dnodes , 1+find(ef_normal==0) );
    end
    
    %-Actual Loop over 2D ips:
    Nf = 0;Af=0;
    for ip2D = 1:n_ips2D
        
        %Compute 2D Interp Data to get detJ_ip:
        [ ~ , ~ , detJ_ip2D , ~ , ~ ] = LagrangeInterp2D_C19N( ef_2Dcoords , Cips2D(ip2D,:) );
        if detJ_ip2D <= 0
            warning('det 2D negative')
        end
        %Evaluate and accumulate elemental H matrix:
        He = He + ( repelem( N3D_ips(ip2D,:) , dofpn ).*repmat( ef_nmat , 1 , npe ) ) * Wips2D(ip2D) * detJ_ip2D;
        
        Af = Af + detJ_ip2D * Wips2D(ip2D);
        Nf = Nf + N3D_ips(ip2D,:) * Wips2D(ip2D) * detJ_ip2D;
    end
    
    % ----------------------------------------------------------------------------- %
    if abs(Af-1/4)>=1e-10
        %disp('Stop here');
    end
    Nf_check = Nf(abs(Nf)>=1e-10);
    if all(Nf(faces_connect(ef,3:end)) == Nf_check)~=1
        disp('Stop here');
    end
    
    %Assembly: 
    H_E( : , : , e ) = H_E( : , : , e ) + He; % ( F I X M E !!! )
    
end

%Trial: vectorize assembly
Hidx_rows = repmat( (1:6)' , 1 , dofpe , n_elem );
Hidx_cols = repmat( permute(dofs_E,[3,2,1]) , 6 , 1 , 1 );

H = sparse( reshape(Hidx_rows,6*dofpe*n_elem,1) ,...
            reshape(Hidx_cols,6*dofpe*n_elem,1) ,...
            reshape(H_E,6*dofpe*n_elem,1) );
%H = sparse( Hidx_rows(:) , Hidx_cols(:) , H_E(:) ); %Slower

return
