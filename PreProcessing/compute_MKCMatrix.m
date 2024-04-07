%clear,clc;load('mesorve3D_masterelem.mat');
function [ H , dofd ] = compute_MKCMatrix( MODEL , quadrature )

MODEL.BoundingFaces.XNormal.Minus.Connectivities = get_face_connectivities( MODEL , MODEL.NodePairs.XNormalFaces(:,1)' );
MODEL.BoundingFaces.XNormal.Plus.Connectivities  = get_face_connectivities( MODEL , MODEL.NodePairs.XNormalFaces(:,2)' );
MODEL.BoundingFaces.XNormal.Minus.Normal = [ -1 ; 0 ; 0 ];
MODEL.BoundingFaces.XNormal.Plus.Normal  = [ 1 ; 0 ; 0 ];

MODEL.BoundingFaces.YNormal.Minus.Connectivities = get_face_connectivities( MODEL , MODEL.NodePairs.YNormalFaces(:,1)' );
MODEL.BoundingFaces.YNormal.Plus.Connectivities  = get_face_connectivities( MODEL , MODEL.NodePairs.YNormalFaces(:,2)' );
MODEL.BoundingFaces.YNormal.Minus.Normal = [ 0 ; -1 ; 0 ];
MODEL.BoundingFaces.YNormal.Plus.Normal  = [ 0 ; 1 ; 0 ];

MODEL.BoundingFaces.ZNormal.Minus.Connectivities = get_face_connectivities( MODEL , MODEL.NodePairs.ZNormalFaces(:,1)' );
MODEL.BoundingFaces.ZNormal.Plus.Connectivities  = get_face_connectivities( MODEL , MODEL.NodePairs.ZNormalFaces(:,2)' );
MODEL.BoundingFaces.ZNormal.Minus.Normal = [ 0 ; 0 ; -1 ];
MODEL.BoundingFaces.ZNormal.Plus.Normal  = [ 0 ; 0 ; 1 ];

options.quadrature = quadrature;


faces_connectivities = [ 1*ones(size(MODEL.BoundingFaces.XNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.XNormal.Minus.Connectivities ;
                         2*ones(size(MODEL.BoundingFaces.XNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.XNormal.Plus.Connectivities  ;
                         3*ones(size(MODEL.BoundingFaces.YNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.YNormal.Minus.Connectivities ;
                         4*ones(size(MODEL.BoundingFaces.YNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.YNormal.Plus.Connectivities  ;
                         5*ones(size(MODEL.BoundingFaces.ZNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.ZNormal.Minus.Connectivities ;
                         6*ones(size(MODEL.BoundingFaces.ZNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.ZNormal.Plus.Connectivities  ;
                         ];
                     
faces_normals = [ MODEL.BoundingFaces.XNormal.Minus.Normal' ;
                  MODEL.BoundingFaces.XNormal.Plus.Normal'  ;
                  MODEL.BoundingFaces.YNormal.Minus.Normal' ;
                  MODEL.BoundingFaces.YNormal.Plus.Normal'  ;
                  MODEL.BoundingFaces.ZNormal.Minus.Normal' ;
                  MODEL.BoundingFaces.ZNormal.Plus.Normal'  ;
                  ];
              
% faces_connectivities = [ 1*ones(size(MODEL.BoundingFaces.XNormal.Minus.Connectivities,1),1) MODEL.BoundingFaces.XNormal.Minus.Connectivities ];                   
% faces_normals = [ MODEL.BoundingFaces.XNormal.Minus.Normal' ];

%Compute BC matrix:
H = Hexa_MKCmat( MODEL , faces_connectivities , faces_normals , options );

%Get 3 (random) dependent dofs: <---- TODOOOOOOO
% nodesd = [MODEL.BoundingFaces.XNormal.Minus.Connectivities(1,1);MODEL.BoundingFaces.ZNormal.Plus.Connectivities(1,1)];
% dofpn = 3;
% dofd   = dofpn*(repelem( nodesd , dofpn , 1 )-1) + repmat((1:dofpn)',size(nodesd,1),1); %dependent dofs
abs_h_mat = abs(full(H));
dofd = zeros(6,1);
for i=1:6
    [ ~ , dofd(i) ] = max(abs_h_mat(i,:)); 
    abs_hmat(i,dofd(i)) = NaN;
end

end

function FaceConect = get_face_connectivities( MODEL , facenodes )

%get model data:
npe    = size(MODEL.Conectivity(:,3:end),2);
nelem  = size(MODEL.Conectivity,1);
if npe == 27 ; npf = 9; end %Nodes per face (adhoc for 27n hexahedra)

%reshape connectivities in column fashion:
resh_elem_per_node  = reshape( (1:nelem)' .* ones(1,npe) , nelem*npe , 1 );
resh_node_per_elem  = reshape( (1:npe) .* ones(nelem,1) , nelem*npe , 1 );
resh_connectivities = reshape( MODEL.Conectivity(:,3:end) , nelem*npe , 1 );

%get elements with at least 1 node on this face:
[ idx_facenodes , ~ ] = find( resh_connectivities == facenodes );
face_elems      = resh_elem_per_node(idx_facenodes);
face_localnodes = resh_node_per_elem(idx_facenodes);

%Loop only over the elements having at least one node on this face:
unique_elems_face = unique( face_elems );
FaceConect        = zeros(length(unique_elems_face), npf+1);

for i=1:length(unique_elems_face)
    
    %get element:
    FaceConect(i,1) = unique_elems_face(i); 
    
    %get local nodes for this element:
    lnodes = sort(face_localnodes( face_elems == unique_elems_face(i) )); %local nodes(sorted)
    if length(lnodes)==npf %if all the nodes in one face are in the boundary, store:
        FaceConect(i,2:end) = lnodes;
    else %leave the row in zeros
        FaceConect(i,:) = 0;
    end
    
end
%Remove (if there are) elements without an entire face on this boundary:
FaceConect( all(FaceConect,2)==0 , : ) = [];

end

