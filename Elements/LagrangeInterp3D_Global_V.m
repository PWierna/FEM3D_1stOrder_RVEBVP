function [ N , detJ , J , invJ , dNGlobal ] = LagrangeInterp3D_Global_V( globalips_conect , InterpData , MODEL )


%% 1ST COMPUTATIONS ============================================================== %%
n_gips   = size(globalips_conect,1); %Number of global integration points
con_gips = MODEL.Conectivity(globalips_conect(:,1),3:end);  %Global Node numbers  

% Corresponding nodal coordinates (Cartesian) for each ip:
% Arrangement:
%   Row     : Global IP N°
%   Col     : [ CoordX CoordY CoordZ ]
%   3rdDim  : Node number (1 to 8)

npe = size(MODEL.Conectivity(:,3:end),2);
n_coords = zeros(n_gips,3,npe); 
n_coords(:,1,:) = reshape(MODEL.Coordinates(con_gips,2),n_gips,1,npe);  %X Coords
n_coords(:,2,:) = reshape(MODEL.Coordinates(con_gips,3),n_gips,1,npe);  %Y Coords
n_coords(:,3,:) = reshape(MODEL.Coordinates(con_gips,4),n_gips,1,npe);  %Z Coords
 

%% PARSE INTERPOLATION DATA ====================================================== %%
% Rearrangement as column vectors: 
%   Row     : Global IP Number
%   3rd Dim : ith node (1 to 8)
N        = permute( InterpData(:,1,globalips_conect(:,2)) , [3,2,1] );
dNdXi    = permute( InterpData(:,2,globalips_conect(:,2)) , [3,2,1] );
dNdEta   = permute( InterpData(:,3,globalips_conect(:,2)) , [3,2,1] );
dNdDseta = permute( InterpData(:,4,globalips_conect(:,2)) , [3,2,1] );


%% COMPUTE JACOBIAN'S COEFFICIENTS FOR EACH IP =================================== %%
J = zeros( n_gips , 9 );
% J(:,1) = sum( dNdXi    .* n_coords(:,:,1) , 3 ); %J11
% J(:,2) = sum( dNdEta   .* n_coords(:,:,1) , 3 ); %J21
% J(:,3) = sum( dNdDseta .* n_coords(:,:,1) , 3 ); %J31
% J(:,4) = sum( dNdXi    .* n_coords(:,:,2) , 3 ); %J12
% J(:,5) = sum( dNdEta   .* n_coords(:,:,2) , 3 ); %J22
% J(:,6) = sum( dNdDseta .* n_coords(:,:,2) , 3 ); %J32
% J(:,7) = sum( dNdXi    .* n_coords(:,:,3) , 3 ); %J13
% J(:,8) = sum( dNdEta   .* n_coords(:,:,3) , 3 ); %J23
% J(:,9) = sum( dNdDseta .* n_coords(:,:,3) , 3 ); %J33

J(:,[1,4,7]) = sum( dNdXi    .* n_coords , 3 ); % [ J11 , J12 , J13 ]
J(:,[2,5,8]) = sum( dNdEta   .* n_coords , 3 ); % [ J21 , J22 , J23 ]
J(:,[3,6,9]) = sum( dNdDseta .* n_coords , 3 ); % [ J31 , J32 , J33 ]


%% COMPUTE JACOBIAN'S DETERMINANT FOR EACH IP ==================================== %%
detJ = J(:,1).*J(:,5).*J(:,9) - J(:,1).*J(:,8).*J(:,6) - J(:,4).*J(:,2).*J(:,9) +...
       J(:,4).*J(:,8).*J(:,3) + J(:,7).*J(:,2).*J(:,6) - J(:,7).*J(:,5).*J(:,3); 

if any(detJ<=0.0)
    warning('Found non-positive Jacobian determinant');
    %disp(detJ_gips(detJ_gips<=0.0));
end
   
if nargout>3
    %% COMPUTE JACOBIAN INVERSE'S COEFFICIENTS FOR EACH IP =========================== %%
    invJ = zeros( n_gips , 9 );
    invJ(:,1) = (1./detJ) .* ( J(:,5).*J(:,9) - J(:,8).*J(:,6) ); %invJ11 
    invJ(:,2) = (1./detJ) .* ( J(:,8).*J(:,3) - J(:,2).*J(:,9) ); %invJ21 
    invJ(:,3) = (1./detJ) .* ( J(:,2).*J(:,6) - J(:,5).*J(:,3) ); %invJ31  
    invJ(:,4) = (1./detJ) .* ( J(:,7).*J(:,6) - J(:,4).*J(:,9) ); %invJ12  
    invJ(:,5) = (1./detJ) .* ( J(:,1).*J(:,9) - J(:,7).*J(:,3) ); %invJ22  
    invJ(:,6) = (1./detJ) .* ( J(:,4).*J(:,3) - J(:,1).*J(:,6) ); %invJ32  
    invJ(:,7) = (1./detJ) .* ( J(:,4).*J(:,8) - J(:,7).*J(:,5) ); %invJ13  
    invJ(:,8) = (1./detJ) .* ( J(:,7).*J(:,2) - J(:,1).*J(:,8) ); %invJ23  
    invJ(:,9) = (1./detJ) .* ( J(:,1).*J(:,5) - J(:,4).*J(:,2) ); %invJ33  
    
    if nargout>4
        %% COMPUTE CARTESIAN DERIVATES =================================================== %%
        %Recall that:
        % dNdX = sum( [ dNdXi , dNdEta , dNdDseta ] .* [ invJ11 invJ12 invJ13 ] , 2 )
        % dNdY = sum( [ dNdXi , dNdEta , dNdDseta ] .* [ invJ21 invJ22 invJ23 ] , 2 )
        % dNdZ = sum( [ dNdXi , dNdEta , dNdDseta ] .* [ invJ31 invJ32 invJ33 ] , 2 )
        
        dNdX = sum( [ dNdXi , dNdEta , dNdDseta ] .* invJ(:,[1,4,7]) , 2 );
        dNdY = sum( [ dNdXi , dNdEta , dNdDseta ] .* invJ(:,[2,5,8]) , 2 );
        dNdZ = sum( [ dNdXi , dNdEta , dNdDseta ] .* invJ(:,[3,6,9]) , 2 );
        
        dNGlobal = [ dNdX , dNdY , dNdZ ];
    end
    
end

end    