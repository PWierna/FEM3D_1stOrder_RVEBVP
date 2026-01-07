function paired_nodes_idx = pairupFaceNodes_nodebynode( ncoords_face1 , ncoords_face2 )

% --------------------------------------------------------------------------------- %
% Inputs:
%  - ncoords_face1 = [ x1 x2 ] (nx2) matrix with the 1-2 coordinates for the n nodes 
%                              lying on the 'master' face.
%  - ncoords_face2 = [ x1 x2 ] (nx2) matrix with the 1-2 coordinates for the n nodes 
%                              lying on the 'slave' face.
% Output:
%  - paired_nodes_idx = (nx1) vector of integers, where the i-th value is the index
%                        of the node in the slave face that is paired to
%                        the i-th node of the master face. That is, the
%                        i-th pair of nodes is [ face1(i) , face2(paired_nodes_idx(i)) ]
% --------------------------------------------------------------------------------- %

%Hard settings - tolerances:
rel_tol = 1e-5;
tol     = rel_tol*max(max(abs(ncoords_face1)));

% Initialize:
paired_nodes_idx = zeros( size(ncoords_face1,1) , 1 );

% Loop over nodes of the master face:
for n = 1:size(ncoords_face1,1)
    
    % Compute in-plane distance from all nodes of the slave face to node n of the
    % master face:
    dist2n  = vecnorm( ncoords_face2 - ncoords_face1(n,:) , 2 , 2 ); 
    
    % Retrieve pair as the node of the slave face with minimum distance:
    [ ~ , paired_nodes_idx(n) ] = min( dist2n );
    
end

% Check both coordinates prior to retrieving pairs:
check1 = ( abs( ncoords_face1(:,1) - ncoords_face2(paired_nodes_idx,1) )<=tol );
check2 = ( abs( ncoords_face1(:,2) - ncoords_face2(paired_nodes_idx,2) )<=tol );

% Report and return:
if all(check1)*all(check2) ~=1
    warning('Node-by-node pairing algorithm was not able to detect some node pairs within the specified tolerance.');
else
    disp('Node-by-node pairing algorithm succeeded.');
end


end