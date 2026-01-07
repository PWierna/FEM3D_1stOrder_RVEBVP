function [ L, dofd ] = build_RVE3D_BCmatrix( MODEL, options )

% "options" structure:
% options.kinadm_conditions_type  = "Periodic" / "Linear" / "Minimal" / "Taylor"
% options.quadrature = "Gauss 3x3x3";

%Hard Settings:
faces_ids = [ "XNormalFaces" ; "YNormalFaces" ; "ZNormalFaces" ];
ndofs  = 3*size(MODEL.Coordinates,1); %(Hard setting: 3 dofpn)

%Get/merge all boundary nodes:
bound_nodes = unique( [ reshape( MODEL.NodePairs.( faces_ids(1) ), [], 1);
                        reshape( MODEL.NodePairs.( faces_ids(2) ), [], 1);
                        reshape( MODEL.NodePairs.( faces_ids(3) ), [], 1) ] );


%% COMPUTE LINEAR BCs MATRIX ===================================================== %%

switch options.kinadm_conditions_type
    case "Taylor"
        %Fix all nodes:
        [ Ibc, dofd ] = buildTaylorBCsmatrix( MODEL );        
        
    case "Linear"
        %Fix all boundary nodes:            
        [ Ibc, dofd ] = buildLinearBCsmatrix( bound_nodes, ndofs );        
        
    case "Periodic"
        
        %Compute periodic bcs:
        %[ Ipbc, dofd_pbc ] = buildPeriodicBCsmatrix( MODEL.NodePairs, ndofs, MODEL.Coordinates(:,2:end) );
        [ Ipbc , dofd_pbc ] = buildPeriodicBCsmatrix( MODEL.NodePairs, ndofs );
        
        % Three additional restrictions are neccessary to prevent rigid body motion, 
        % thus we fix one interior node (randomly-chosen as the first node in the set
        % of interior nodes):
        interior_nodes = setdiff( MODEL.Coordinates(:,1), bound_nodes );
        
        [ Ibc_essential,...
          dofd_essential ] = buildEssentialBCsmatrix( interior_nodes(1),... %x-fixed node
                                                      interior_nodes(1),... %y-fixed node
                                                      interior_nodes(1),... %z-fixed node
                                                      ndofs );
        
        Ibc  = [ Ipbc ; Ibc_essential ];
        dofd = [ dofd_pbc ; dofd_essential ];
        
        
    case "Minimal"
        %Compute minimal kinematic constraints matrix:
        [ Imkc , dofd_mkc ] = compute_MKCMatrix( MODEL , options.quadrature ); 
        
        % Five additional (independent) restrictions are neccessary to prevent rigid 
        % body motion and rigid body rotation. Thus we fix 2 interior nodes 
        % accordingly:
        interior_nodes = setdiff( MODEL.Coordinates(:,1), bound_nodes );

        [ Ibc_essential,...
          dofd_essential ] = buildEssentialBCsmatrix( interior_nodes(1:2) ,... %x-fixed nodes
                                                      interior_nodes(1:2) ,... %y-fixed nodes
                                                      interior_nodes(1) ,...   %z-fixed node
                                                      ndofs );
                                                                  
        Ibc  = [ Imkc ; Ibc_essential ];
        dofd = [ dofd_mkc ; dofd_essential ];
        
end


% TO-DO: Sort rows to reduce bandwidth 
% [ ~ , sortidx ] = sort(dofd);
% dofd = dofd(sortidx);
% doff = (1:ndofs)'; doff(dofd)=[];
% Ibc  = Ibc(sortidx,:);
% faces_nodes_all=unique(faces_nodes_all); 
% [ Ibc , dofd ] = buildLinearBCsmatrix( faces_nodes_all ,  ndofs );

% Compute Dofs vectors:
dofd = sort(dofd); 
doff = (1:ndofs)'; doff(dofd)=[];

% Compute partition and retrieve Linear BCs Matrix:
Id = Ibc(:,dofd); If = Ibc(:,doff);
L  = - Id \ If;

end


%% ROUTINES TO ASSEMBLE BCs MATRICES ============================================= %%

%ESSENTIAL BCs MATRIX
function [ BCmat, dofd ] = buildEssentialBCsmatrix( fixed_x, fixed_y, fixed_z, ndofs )

nc     = length([fixed_x; fixed_y; fixed_z]); %number of constraints 
idxrow = (1:nc)';                             %row indexes = constraint number 
idxcol = [ 3*( fixed_x - 1 ) + 1 ;
           3*( fixed_y - 1 ) + 2 ;
           3*( fixed_z - 1 ) + 3 ];           %col indexes = ndof 
       

% Assemble essential BC matrix
BCmat = sparse( idxrow , idxcol , ones(nc,1) , nc , ndofs );

% Dependent dofs
dofd = idxcol;       

end

% LINEAR BCs MATRIX
function [ BCsmat, dofd ] = buildLinearBCsmatrix( faces_nodes, ndofs )

%Assert there are no repeated nodes:
faces_nodes = unique(faces_nodes);

%Auxiliar computations:
nc      = 3*length(faces_nodes);                 %number of constraints
idxrow  = (1:nc)';                               %row indexes = constraint number
idxcol  = reshape(3*(faces_nodes-1)+(1:3),nc,1); %col indexes = ndof associated w/each constraint

%Linear BC matrix:
BCsmat = sparse( idxrow, idxcol, ones(nc,1), nc, ndofs );

%Dependent dofs:
dofd = idxcol;

end

%TAYLOR BCs MATRIX
function [ BCmat , dofd ] = buildTaylorBCsmatrix( MODEL )

ndofs  = 3*size( MODEL.Coordinates, 1 ); %Get total number of dofs

%Taylor BC matrix:
BCmat = sparse(eye(ndofs));

%Dependent dofs:
dofd = (1:ndofs)';

end

% PERIODIC BCs MATRIX
function [ BCmat, dofd ] = buildPeriodicBCsmatrix( NodePairs, ndofs )

% Helper function to extract nodes from a given set of nodes
function diffnodes = removenodepairs( mainset, nodes2extract )
    
%Indexes of the mainset rows whose (all) elements are not in the set to be
%extracted:
idx2keep = all( ~ismember(mainset,nodes2extract), 2 );

%Return only the rows not having nodes that are in the set to extract:
diffnodes = mainset(idx2keep,:);

end

%Hard setting (3D Problem)
dofpn = 3;

% Get nodes corresponding to the eight RVE-vertices:
vertexnodes.xminyminzmin = intersect( intersect(NodePairs.XNormalFaces(:,1), NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xminyminzmax = intersect( intersect(NodePairs.XNormalFaces(:,1), NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.xminymaxzmin = intersect( intersect(NodePairs.XNormalFaces(:,1), NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xminymaxzmax = intersect( intersect(NodePairs.XNormalFaces(:,1), NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.xmaxyminzmin = intersect( intersect(NodePairs.XNormalFaces(:,2), NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xmaxyminzmax = intersect( intersect(NodePairs.XNormalFaces(:,2), NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.xmaxymaxzmin = intersect( intersect(NodePairs.XNormalFaces(:,2), NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xmaxymaxzmax = intersect( intersect(NodePairs.XNormalFaces(:,2), NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.all = [ vertexnodes.xminyminzmin; vertexnodes.xminyminzmax;
                    vertexnodes.xminymaxzmin; vertexnodes.xminymaxzmax;
                    vertexnodes.xmaxyminzmin; vertexnodes.xmaxyminzmax;
                    vertexnodes.xmaxymaxzmin; vertexnodes.xmaxymaxzmax ];

% Get edges' "internal" nodes, i.e. nodes lying at the RVE's edges excluding its 
% (vertex) boundaries:           
edgesintnodes.xminymin = removenodepairs( intersect(NodePairs.XNormalFaces(:,1), NodePairs.YNormalFaces(:,1)), vertexnodes.all );
edgesintnodes.xminymax = removenodepairs( intersect(NodePairs.XNormalFaces(:,1), NodePairs.YNormalFaces(:,2)), vertexnodes.all );
edgesintnodes.xminzmin = removenodepairs( intersect(NodePairs.XNormalFaces(:,1), NodePairs.ZNormalFaces(:,1)), vertexnodes.all );
edgesintnodes.xminzmax = removenodepairs( intersect(NodePairs.XNormalFaces(:,1), NodePairs.ZNormalFaces(:,2)), vertexnodes.all );
edgesintnodes.xmaxymin = removenodepairs( intersect(NodePairs.XNormalFaces(:,2), NodePairs.YNormalFaces(:,1)), vertexnodes.all );
edgesintnodes.xmaxymax = removenodepairs( intersect(NodePairs.XNormalFaces(:,2), NodePairs.YNormalFaces(:,2)), vertexnodes.all );
edgesintnodes.xmaxzmin = removenodepairs( intersect(NodePairs.XNormalFaces(:,2), NodePairs.ZNormalFaces(:,1)), vertexnodes.all );
edgesintnodes.xmaxzmax = removenodepairs( intersect(NodePairs.XNormalFaces(:,2), NodePairs.ZNormalFaces(:,2)), vertexnodes.all );
edgesintnodes.yminzmin = removenodepairs( intersect(NodePairs.YNormalFaces(:,1), NodePairs.ZNormalFaces(:,1)), vertexnodes.all );
edgesintnodes.yminzmax = removenodepairs( intersect(NodePairs.YNormalFaces(:,1), NodePairs.ZNormalFaces(:,2)), vertexnodes.all );
edgesintnodes.ymaxzmin = removenodepairs( intersect(NodePairs.YNormalFaces(:,2), NodePairs.ZNormalFaces(:,1)), vertexnodes.all );
edgesintnodes.ymaxzmax = removenodepairs( intersect(NodePairs.YNormalFaces(:,2), NodePairs.ZNormalFaces(:,2)), vertexnodes.all );
edgesintnodes.all = [ edgesintnodes.xminymin; edgesintnodes.xminymax; 
                      edgesintnodes.xminzmin; edgesintnodes.xminzmax;
                      edgesintnodes.xmaxymin; edgesintnodes.xmaxymax;
                      edgesintnodes.xmaxzmin; edgesintnodes.xmaxzmax;
                      edgesintnodes.yminzmin; edgesintnodes.yminzmax;
                      edgesintnodes.ymaxzmin; edgesintnodes.ymaxzmax ];

% Get "interior" bounding faces nodes, i.e. nodes lying at the RVE's faces excluding 
%its (edges and vertices) boundaries:
bndfacesintnodes.xnormal = [ removenodepairs(NodePairs.XNormalFaces(:,1), [edgesintnodes.all;vertexnodes.all]) ,... %face xmin
                             removenodepairs(NodePairs.XNormalFaces(:,2), [edgesintnodes.all;vertexnodes.all]) ];   %face xmax
bndfacesintnodes.ynormal = [ removenodepairs(NodePairs.YNormalFaces(:,1), [edgesintnodes.all;vertexnodes.all]) ,... %face ymin
                             removenodepairs(NodePairs.YNormalFaces(:,2), [edgesintnodes.all;vertexnodes.all]) ];   %face ymax
bndfacesintnodes.znormal = [ removenodepairs(NodePairs.ZNormalFaces(:,1), [edgesintnodes.all;vertexnodes.all]) ,... %face zmin
                             removenodepairs(NodePairs.ZNormalFaces(:,2), [edgesintnodes.all;vertexnodes.all]) ];   %face zmax
                            
% Merge all "interior" bounding faces pairs of nodes 
% ([normal<0=dependent , normal>0=free])
bndfacesintnodes.all = [ bndfacesintnodes.xnormal ; 
                         bndfacesintnodes.ynormal ;
                         bndfacesintnodes.znormal ];

% Assemble node pairs for "interior" edges nodes ([Dependent, Free]). Note
% that there are three edges, [xmin,ymin,z], [xmin,y,zmin], and [x,ymin,zmin] whose 
% nodes will be taken as master (free) nodes of the nodes in the remaining 9 
% (dependent/slave) edges of the prismatic RVE, as per the following convention:
%
%   Pair       Slave --------------> Master
%    1     [xmax,ymin,z] ------>  [xmin,ymin,z]  (Z-Oriented)
%    2     [xmin,ymax,z] ------>  [xmin,ymin,z]  (Z-Oriented)
%    3     [xmax,ymax,z] ------>  [xmin,ymin,z]  (Z-Oriented)
%    4     [xmin,y,zmax] ------>  [xmin,y,zmin]  (Y-Oriented)
%    5     [xmax,y,zmin] ------>  [xmin,y,zmin]  (Y-Oriented)
%    6     [xmax,y,zmax] ------>  [xmin,y,zmin]  (Y-Oriented)
%    7     [x,ymin,zmax] ------>  [x,ymin,zmin]  (X-Oriented)
%    8     [x,ymax,zmin] ------>  [x,ymin,zmin]  (X-Oriented)
%    9     [x,ymax,zmax] ------>  [x,ymin,zmin]  (X-Oriented)

% Pair 1 [ (xmax,ymin,z), (xmin,ymin,z) ]:
edgespair1 = zeros( length(edgesintnodes.xminymin) , 2 );
edgespair1(:,2) = edgesintnodes.xminymin;
for i = 1:length( edgesintnodes.xminymin )
    % Look for pairs using the opposite X-Normal face(Xmax)
    edgespair1(i,1) = NodePairs.XNormalFaces( NodePairs.XNormalFaces(:,1) == edgesintnodes.xminymin(i) , 2 );
end

% Pair 2 [ (xmin,ymax,z), (xmin,ymin,z) ]:
edgespair2 = zeros( length(edgesintnodes.xminymin) , 2 );
edgespair2(:,2) = edgesintnodes.xminymin;
for i = 1:length( edgesintnodes.xminymin )
    % Look for pairs using the opposite Y-Normal face(Ymax)
    edgespair2(i,1) = NodePairs.YNormalFaces( NodePairs.YNormalFaces(:,1) == edgesintnodes.xminymin(i) , 2 );
end

% Edges Pair 3 [ (xmax,ymax,z), (xmin,ymin,z) ]:
edgespair3      = zeros( length(edgesintnodes.xminymin) , 2 );
edgespair3(:,2) = edgesintnodes.xminymin;
for i = 1:length( edgesintnodes.xminymin )
    % Look for pairs using pair 1 and the opposite Y-Normal face (Ymax)
    edgespair3(i,1) = NodePairs.YNormalFaces( NodePairs.YNormalFaces(:,1) == edgespair1(i,1) , 2 );
end

% Edges Pair 4 [ (xmin,y,zmax), (xmin,y,zmin) ]:
edgespair4      = zeros( length(edgesintnodes.xminzmin) , 2 );
edgespair4(:,2) = edgesintnodes.xminzmin;
for i = 1:length(edgesintnodes.xminzmin)
    % Look for pairs in the Z-Normal opposite face (ZMax):
    edgespair4(i,1) = NodePairs.ZNormalFaces( NodePairs.ZNormalFaces(:,1) == edgesintnodes.xminzmin(i) , 2 );
end

% Edges Pair 5 [ (xmax,y,zmin), (xmin,y,zmin) ]:
edgespair5      = zeros( length(edgesintnodes.xminzmin) , 2 );
edgespair5(:,2) = edgesintnodes.xminzmin;
for i = 1:length(edgesintnodes.xminzmin)
    % Look for pairs in the X-Normal opposite face (XMax):
    edgespair5(i,1) = NodePairs.XNormalFaces( NodePairs.XNormalFaces(:,1) == edgesintnodes.xminzmin(i) , 2 );
end

% Edges Pair 6 [ (xmax,y,zmax), (xmin,y,zmin) ]:
edgespair6      = zeros( length(edgesintnodes.xminzmin) , 2 );
edgespair6(:,2) = edgesintnodes.xminzmin;
for i = 1:length(edgesintnodes.xminzmin)
    % Look for pairs using pair 5 and the opposite Z-Normal face (Zmax):
    edgespair6(i,1) = NodePairs.ZNormalFaces( NodePairs.ZNormalFaces(:,1)==edgespair5(i,1) , 2 );
end

% Edges Pair 7 [ (x,ymin,zmax), (x,ymin,zmin) ]:
edgespair7      = zeros( length(edgesintnodes.yminzmin) , 2 );
edgespair7(:,2) = edgesintnodes.yminzmin;
for i = 1:length(edgesintnodes.yminzmin)
    % Look for pairs in the Z-Normal opposite face (ZMax):
    edgespair7(i,1) = NodePairs.ZNormalFaces( NodePairs.ZNormalFaces(:,1)==edgesintnodes.yminzmin(i) , 2 );
end

% Edges Pair 8 [ (x,ymax,zmin), (x,ymin,zmin) ]:
edgespair8      = zeros( length(edgesintnodes.yminzmin) , 2 );
edgespair8(:,2) = edgesintnodes.yminzmin;
for i = 1:length(edgesintnodes.yminzmin)
    % Look for pairs in the Y-Normal opposite face (YMax):
    edgespair8(i,1) = NodePairs.YNormalFaces( NodePairs.YNormalFaces(:,1)==edgesintnodes.yminzmin(i) , 2 );
end

% Edges Pair 9 [ (x,ymax,zmax), (x,ymin,zmin) ]:
edgespair9      = zeros( length(edgesintnodes.yminzmin) , 2 );
edgespair9(:,2) = edgesintnodes.yminzmin;
for i = 1:length(edgesintnodes.yminzmin)
    % Look for pairs using pair 7 and the opposite Y-Normal face (Ymax):
    edgespair9(i,1) = NodePairs.YNormalFaces( NodePairs.YNormalFaces(:,1)==edgespair7(i,1) , 2 );
end

% Concatentate edge's node pairs
npairs_intedges = [ edgespair1; edgespair2; edgespair3; edgespair4;
                    edgespair5; edgespair6; edgespair7; edgespair8; edgespair9 ];

% Assemble node pairs for vertex nodes ( [Dependent, Free] ):
npairs_vertex = [ vertexnodes.xminyminzmax, vertexnodes.xminyminzmin;
                  vertexnodes.xminymaxzmin, vertexnodes.xminyminzmin;
                  vertexnodes.xminymaxzmax, vertexnodes.xminyminzmin;
                  vertexnodes.xmaxyminzmin, vertexnodes.xminyminzmin;
                  vertexnodes.xmaxyminzmax, vertexnodes.xminyminzmin;
                  vertexnodes.xmaxymaxzmin, vertexnodes.xminyminzmin;
                  vertexnodes.xmaxymaxzmax, vertexnodes.xminyminzmin ];

%Concatenate all pair of nodes( [ dependent, free ] ):                                  
npairs = [ bndfacesintnodes.all ; npairs_intedges ; npairs_vertex ]; 

%Compute dofs for all the pairs of nodes:
nconstr = dofpn * size(npairs,1); %number of constraints
dofd    = dofpn*(repelem( npairs(:,1) , dofpn , 1 )-1) + repmat((1:dofpn)',size(npairs,1),1); %dependent dofs
doff    = dofpn*(repelem( npairs(:,2) , dofpn , 1 )-1) + repmat((1:dofpn)',size(npairs,1),1); %free dofs

%Assemble Periodic BC matrix:
idxrow  = [ (1:nconstr)' ; (1:nconstr)' ];        
idxcol  = [ dofd ; doff ];              
cvals   = [ ones(nconstr,1) ; -ones(nconstr,1) ]; %constraints values

BCmat = sparse( idxrow , idxcol , cvals  , nconstr , ndofs );

end

