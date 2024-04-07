function [ L , dofd ] = build_RVE3D_BCmatrix( MODEL , options )

%Hard Settings:
faces_ids = [ "XNormalFaces" ; "YNormalFaces" ; "ZNormalFaces" ];
ndofs  = 3*size(MODEL.Coordinates,1); %(Hard setting: 3 dofpn)

% "options" structure:
% options.kinadm_conditions_type  = "Periodic" / "Linear" / "Minimal" / "Taylor"
% options.quadrature = "Gauss 3x3x3";

%Get/merge all boundary nodes:
bound_nodes = unique( [ reshape( MODEL.NodePairs.( faces_ids(1) ) , [] , 1 );
                        reshape( MODEL.NodePairs.( faces_ids(2) ) , [] , 1 );
                        reshape( MODEL.NodePairs.( faces_ids(3) ) , [] , 1 ) ] );


%% COMPUTE BC MATRIX ============================================================= %%


switch options.kinadm_conditions_type
    case "Taylor"
        %Fix all nodes:
        [ Ibc , dofd ] = build_TaylorBCmatrix( MODEL );
        
        
    case "Linear"
        %Fix all boundary nodes:            
        [ Ibc , dofd ] = build_LinearBCmatrix( bound_nodes , ndofs );
        
        
    case "Periodic"
        %Compute periodic bcs:
        [ Ipbc , dofd_pbc ] = build_PeriodicBCmatrix( MODEL.NodePairs , ndofs );
        
        %3 additional restrictions are neccessary to prevent rigid body motion, thus we
        %fix one interior node:
        interior_nodes = setdiff( MODEL.Coordinates(:,1) , bound_nodes );
        
        [ Ibc_essential , dofd_essential ] = build_EssentialBCmatrix( interior_nodes(1),... %x-fixed node
                                                                      interior_nodes(1),... %y-fixed node
                                                                      interior_nodes(1),... %z-fixed node
                                                                      ndofs );
        
        Ibc  = [ Ipbc ; Ibc_essential ];
        dofd = [ dofd_pbc ; dofd_essential ];
        
        
    case "Minimal"
        %Compute minimal kinematic constraints matrix:
        [ Imkc , dofd_mkc ] = compute_MKCMatrix( MODEL , options.quadrature ); 
        
        %5 additional (independent) restrictions are neccessary to prevent rigid body
        %motion and rigid body rotation. Thus we -accordingly- fix 2 interior nodes:
        interior_nodes = setdiff( MODEL.Coordinates(:,1) , bound_nodes );

        [ Ibc_essential , dofd_essential ] = build_EssentialBCmatrix( interior_nodes(1:2) ,... %x-fixed nodes
                                                                      interior_nodes(1:2) ,... %y-fixed nodes
                                                                      interior_nodes(1) ,...   %z-fixed node
                                                                      ndofs );
                                                                  
        Ibc  = [ Imkc ; Ibc_essential ];
        dofd = [ dofd_mkc ; dofd_essential ];
        
end


% %Sort rows to reduce bandwidth: TO-DO!!!!!
% [ ~ , sortidx ] = sort(dofd);
% dofd = dofd(sortidx);
% doff = (1:ndofs)'; doff(dofd)=[];
% Ibc  = Ibc(sortidx,:);
% faces_nodes_all=unique(faces_nodes_all); 
% [ Ibc , dofd ] = build_LinearBCmatrix( faces_nodes_all ,  ndofs );

%Compute Dofs vectors:
dofd = sort(dofd); 
doff = (1:ndofs)'; doff(dofd)=[];

%Partition and compute Linear BC Matrix:
Id = Ibc(:,dofd); If = Ibc(:,doff);
L  = - Id \ If;

end


%% BC MATRICES FUNCTIONS ========================================================= %%

%LINEAR BC MATRIX:
function [ BCmat , dofd ] = build_LinearBCmatrix( faces_nodes , ndofs )

%Assert there are no repeated nodes:
faces_nodes = unique(faces_nodes);

%Auxiliar computations:
nc      = 3*length(faces_nodes);                    %number of constraints
idxrow  = (1:nc)';                                  %row indexes = constraint number
idxcol  = reshape(3*(faces_nodes-1)+(1:3),nc,1);    %col indexes = ndof associated w/each constraint

%Linear BC matrix:
BCmat = sparse( idxrow , idxcol , ones(nc,1) , nc , ndofs );

%Dependent dofs:
dofd = idxcol;

end

%TAYLOR BC MATRIX:
function [ BCmat , dofd ] = build_TaylorBCmatrix( MODEL )

ndofs  = 3*size( MODEL.Coordinates , 1 );   %Get total number of dofs
% nc     = ndofs;                             %number of constraints (= all dofs)
% idxrow = (1:nc)';                           %row indexes = constraint number = node number
% idxcol = (1:ndofs)';                        %col indexes = ndof associated w/each constrains
% 
% %Taylor BC matrix:
% BCmat = sparse( idxrow , idxcol , ones(nc,1) , nc , ndofs );
% 
% %Dependent dofs:
% dofd = idxcol;

%Taylor BC matrix:
BCmat = sparse(eye(ndofs));

%Dependent dofs:
dofd = (1:ndofs)';

end

%PERIODIC BC MATRIX:
function [ BCmat , dofd ] = build_PeriodicBCmatrix( NodePairs, ndofs )

function diffnodes = removenodepairs( mainset , nodes2extract )

%Indexes of the mainset rows whose (all) elements are not in the set to be
%extracted:
idx2keep = all( ~ismember(mainset,nodes2extract) , 2 );

%Return only the rows not having nodes that are in the set to extract:
diffnodes = mainset(idx2keep,:);

end


%Hard setting:
dofpn = 3;

%Get vertex nodes:
vertexnodes.xminyminzmin = intersect( intersect(NodePairs.XNormalFaces(:,1),NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xminyminzmax = intersect( intersect(NodePairs.XNormalFaces(:,1),NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.xminymaxzmin = intersect( intersect(NodePairs.XNormalFaces(:,1),NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xminymaxzmax = intersect( intersect(NodePairs.XNormalFaces(:,1),NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.xmaxyminzmin = intersect( intersect(NodePairs.XNormalFaces(:,2),NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xmaxyminzmax = intersect( intersect(NodePairs.XNormalFaces(:,2),NodePairs.YNormalFaces(:,1)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.xmaxymaxzmin = intersect( intersect(NodePairs.XNormalFaces(:,2),NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,1) );
vertexnodes.xmaxymaxzmax = intersect( intersect(NodePairs.XNormalFaces(:,2),NodePairs.YNormalFaces(:,2)), NodePairs.ZNormalFaces(:,2) );
vertexnodes.all = [ vertexnodes.xminyminzmin;vertexnodes.xminyminzmax;vertexnodes.xminymaxzmin;vertexnodes.xminymaxzmax ;
                    vertexnodes.xmaxyminzmin;vertexnodes.xmaxyminzmax;vertexnodes.xmaxymaxzmin;vertexnodes.xmaxymaxzmax ];

%Get edges "internal" nodes:           
edgesintnodes.xminymin = removenodepairs( intersect(NodePairs.XNormalFaces(:,1),NodePairs.YNormalFaces(:,1)) , vertexnodes.all );
edgesintnodes.xminymax = removenodepairs( intersect(NodePairs.XNormalFaces(:,1),NodePairs.YNormalFaces(:,2)) , vertexnodes.all );
edgesintnodes.xminzmin = removenodepairs( intersect(NodePairs.XNormalFaces(:,1),NodePairs.ZNormalFaces(:,1)) , vertexnodes.all );
edgesintnodes.xminzmax = removenodepairs( intersect(NodePairs.XNormalFaces(:,1),NodePairs.ZNormalFaces(:,2)) , vertexnodes.all );
edgesintnodes.xmaxymin = removenodepairs( intersect(NodePairs.XNormalFaces(:,2),NodePairs.YNormalFaces(:,1)) , vertexnodes.all );
edgesintnodes.xmaxymax = removenodepairs( intersect(NodePairs.XNormalFaces(:,2),NodePairs.YNormalFaces(:,2)) , vertexnodes.all );
edgesintnodes.xmaxzmin = removenodepairs( intersect(NodePairs.XNormalFaces(:,2),NodePairs.ZNormalFaces(:,1)) , vertexnodes.all );
edgesintnodes.xmaxzmax = removenodepairs( intersect(NodePairs.XNormalFaces(:,2),NodePairs.ZNormalFaces(:,2)) , vertexnodes.all );
edgesintnodes.yminzmin = removenodepairs( intersect(NodePairs.YNormalFaces(:,1),NodePairs.ZNormalFaces(:,1)) , vertexnodes.all );
edgesintnodes.yminzmax = removenodepairs( intersect(NodePairs.YNormalFaces(:,1),NodePairs.ZNormalFaces(:,2)) , vertexnodes.all );
edgesintnodes.ymaxzmin = removenodepairs( intersect(NodePairs.YNormalFaces(:,2),NodePairs.ZNormalFaces(:,1)) , vertexnodes.all );
edgesintnodes.ymaxzmax = removenodepairs( intersect(NodePairs.YNormalFaces(:,2),NodePairs.ZNormalFaces(:,2)) , vertexnodes.all );
edgesintnodes.all = [ edgesintnodes.xminymin;edgesintnodes.xminymax;edgesintnodes.xminzmin;edgesintnodes.xminzmax;
                      edgesintnodes.xmaxymin;edgesintnodes.xmaxymax;edgesintnodes.xmaxzmin;edgesintnodes.xmaxzmax;
                      edgesintnodes.yminzmin;edgesintnodes.yminzmax;edgesintnodes.ymaxzmin;edgesintnodes.ymaxzmax ];

%Get "interior" bounding faces nodes:
bndfacesintnodes.xnormal = [ removenodepairs(NodePairs.XNormalFaces(:,1),[edgesintnodes.all;vertexnodes.all]) ,... %face xmin
                             removenodepairs(NodePairs.XNormalFaces(:,2),[edgesintnodes.all;vertexnodes.all]) ];   %face xmax
bndfacesintnodes.ynormal = [ removenodepairs(NodePairs.YNormalFaces(:,1),[edgesintnodes.all;vertexnodes.all]) ,... %face ymin
                             removenodepairs(NodePairs.YNormalFaces(:,2),[edgesintnodes.all;vertexnodes.all]) ];   %face ymax
bndfacesintnodes.znormal = [ removenodepairs(NodePairs.ZNormalFaces(:,1),[edgesintnodes.all;vertexnodes.all]) ,... %face zmin
                             removenodepairs(NodePairs.ZNormalFaces(:,2),[edgesintnodes.all;vertexnodes.all]) ];   %face zmax
                            
%Merge all "interior" bounding faces pairs of nodes ([normal<0=dependent , normal>0=free]):
bndfacesintnodes.all = [ bndfacesintnodes.xnormal ; 
                         bndfacesintnodes.ynormal ;
                         bndfacesintnodes.znormal ];

%Assemble node pairs for "interior" edges nodes ( [Dependent,Free] ):
npairs_intedges = [ edgesintnodes.xmaxymax , edgesintnodes.xminymin ;
                    edgesintnodes.xminymax , edgesintnodes.xminymin ;
                    edgesintnodes.xmaxymin , edgesintnodes.xminymin ;
                    edgesintnodes.xmaxzmax , edgesintnodes.xminzmin ;
                    edgesintnodes.xminzmax , edgesintnodes.xminzmin ;
                    edgesintnodes.xmaxzmin , edgesintnodes.xminzmin ;
                    edgesintnodes.ymaxzmax , edgesintnodes.yminzmin ;
                    edgesintnodes.yminzmax , edgesintnodes.yminzmin ;
                    edgesintnodes.ymaxzmin , edgesintnodes.yminzmin ];

%Assemble node pairs for vertex nodes ( [Dependent, Free] ):
npairs_vertex = [ vertexnodes.xminyminzmax , vertexnodes.xminyminzmin ;
                  vertexnodes.xminymaxzmin , vertexnodes.xminyminzmin ;
                  vertexnodes.xminymaxzmax , vertexnodes.xminyminzmin ;
                  vertexnodes.xmaxyminzmin , vertexnodes.xminyminzmin ;
                  vertexnodes.xmaxyminzmax , vertexnodes.xminyminzmin ;
                  vertexnodes.xmaxymaxzmin , vertexnodes.xminyminzmin ;
                  vertexnodes.xmaxymaxzmax , vertexnodes.xminyminzmin ];

%Concatenate all pair of nodes( [ dependent , free ] ):                                  
npairs = [ bndfacesintnodes.all ; npairs_intedges ; npairs_vertex ]; 

%Compute dofs for all the pairs of nodes:
nconstr = dofpn * size(npairs,1); %number of constraints
dofd    = dofpn*(repelem( npairs(:,1) , dofpn , 1 )-1) + repmat((1:dofpn)',size(npairs,1),1); %dependent dofs
doff    = dofpn*(repelem( npairs(:,2) , dofpn , 1 )-1) + repmat((1:dofpn)',size(npairs,1),1); %free dofs

%Assemble Periodic BC matrix:
idxrow  = [ (1:nconstr)' ; (1:nconstr)' ];        
idxcol  = [ dofd ; doff ];              
cvals   = [ ones(nconstr,1) ; -ones(nconstr,1) ]; %constrains values

BCmat = sparse( idxrow , idxcol , cvals  , nconstr , ndofs );

end

%ESSENTIAL BC MATRIX:
function [ BCmat , dofd ] = build_EssentialBCmatrix( fixed_x , fixed_y , fixed_z , ndofs )

nc     = length( [ fixed_x;fixed_y;fixed_z ] ); %number of constraints 
idxrow = (1:nc)';                               %row indexes = constraint number 
idxcol = [ 3 *( fixed_x - 1 ) + 1 ;
           3 *( fixed_y - 1 ) + 2 ;
           3 *( fixed_z - 1 ) + 3 ];            %col indexes = ndof 
       

%Essential BC matrix:
BCmat = sparse( idxrow , idxcol , ones(nc,1) , nc , ndofs );

%Dependent dofs:
dofd = idxcol;       

end