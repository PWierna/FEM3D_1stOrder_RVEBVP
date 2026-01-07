function MODEL = identifyRVENodePairs( MODEL )
% This script takes a raw RVE-model mesh (usually obtained after reading an
% ANSYS-generated mesh) and gives it the proper structure for the posterior
% RVE-Problem solver. It also generates the data associated to the node
% pairs necessary for the application of Periodic BCs. A prismatic RVE
% model is assumed.

% Retrieve absolute tolerances:
rel_tol = 1e-8;                  % Relative tol (Hard setting)
tol = max(max(MODEL.Coordinates(:,2:4))) * rel_tol;  % Absolute tolerance

% Build node pairs for each pair of faces and allocate in MODEL structure:
MODEL.NodePairs.XNormalFaces = pairupFaceNodes( MODEL.Coordinates(:,2:4) , 1 , tol );
MODEL.NodePairs.YNormalFaces = pairupFaceNodes( MODEL.Coordinates(:,2:4) , 2 , tol );
MODEL.NodePairs.ZNormalFaces = pairupFaceNodes( MODEL.Coordinates(:,2:4) , 3 , tol );




% Uncomment the following lines to display nodes per face
%bounds = [ min(c123)' , max(c123)' ];
%n_c1faces = [ find( abs(c123(:,1)-bounds(1,1)) <= tol ), find( abs(c123(:,1)-bounds(1,2)) <= tol ) ];
%n_c2faces = [ find( abs(c123(:,2)-bounds(2,1)) <= tol ), find( abs(c123(:,2)-bounds(2,2)) <= tol ) ];
%n_c3faces = [ find( abs(c123(:,3)-bounds(3,1)) <= tol ), find( abs(c123(:,3)-bounds(3,2)) <= tol ) ];
%plot3(c123(n_c1faces,1),c123(n_c1faces,2),c123(n_c1faces,3),'or','LineWidth',3,'MarkerFaceColor','r');hold on; grid on; axis equal
%plot3(c123(n_c2faces,1),c123(n_c2faces,2),c123(n_c2faces,3),'ob','LineWidth',3,'MarkerFaceColor','b');
%plot3(c123(n_c3faces,1),c123(n_c3faces,2),c123(n_c3faces,3),'ok','LineWidth',3,'MarkerFaceColor','k');


return

function npairs = pairupFaceNodes( c123 , ci , abstol )

% Get boundary faces along the i-th coordinate 'ci' (a prismatic RVE is assumed):
cibounds = [ min(c123(:,ci)) , max(c123(:,ci)) ];

% Retrieve other coordinates than ci:
otherc = 1:3; otherc(ci) = [];

% Get nodes at the two bounding faces [ ci_min ci_max ] (along the main coordinate):
n_cifaces = [ find( abs(c123(:,ci)-cibounds(1)) <= abstol ), find( abs(c123(:,ci)-cibounds(2)) <= abstol ) ];

% Retrieve paired nodes([ ci_min face , ci_max face ]) by sorting along the other two coordinates:
npairs = pairupFaceNodes_bysorting( n_cifaces(:,2) , n_cifaces(:,1) , c123(:,otherc) );

% Run Check:
checkI  = ( abs(c123(npairs(:,1),otherc(1)) - c123(npairs(:,2),otherc(1))) <= abstol );
checkII = ( abs(c123(npairs(:,1),otherc(2)) - c123(npairs(:,2),otherc(2))) <= abstol );

% If 'pairing by sorting' strategy failed, resort to node-by-node pairing:
if all(checkI)*all(checkII) ~=1
    
    %Report
    fprintf('\nPairing-by-sorting algorithm: Coordinates of node pairs for Faces with Coordinate %d-Normal \nseem not to be within specified tolerance. Switching to node-by-node pairing algorithm...\n',ci);
    
    %Find pairs with node-by-node algortihm
    npairs      = n_cifaces;
    npairs_idx  = pairupFaceNodes_nodebynode( c123(n_cifaces(:,2),otherc) , c123(n_cifaces(:,1),otherc) );
    npairs(:,1) = npairs(npairs_idx,1); %(Reorder to pairup)
    
    %Run check once again
    checkI  = ( abs(c123(npairs(:,1),otherc(1)) - c123(npairs(:,2),otherc(1))) <= abstol );
    checkII = ( abs(c123(npairs(:,1),otherc(2)) - c123(npairs(:,2),otherc(2))) <= abstol );
    if all(checkI)*all(checkII) ~=1
        error('Coordinate '+string(ci)+'-Normal Faces Pairing: Node-by-node pairing algorithm failed. Please ensure that your mesh satisfies periodic constraints.')
    end
    
end