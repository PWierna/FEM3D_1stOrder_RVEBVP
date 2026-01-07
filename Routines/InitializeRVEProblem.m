function [ displ, int_vars ] = InitializeRVEProblem(RVEMODEL, Options)

% --------------------------------------------------------------------------------- % 
% This function initializes both the vector of (unknown) nodal fluctuating 
% displacements and the array of internal variables (at integration points), previous 
% to solving the RVE BVP.
% --------------------------------------------------------------------------------- %

% Initialize nodal displacement vectors
dofpn = 3;	%dofs per node(Hard setting for a 3D problem)
ndofs = dofpn * size(RVEMODEL.Coordinates, 1); %Total N°dofs
displ = zeros(ndofs, 1);


% Get total number of Integration points in the mesh (NOTE: It is assumed that all 
% elements use the same quadrature) and initialize internal variables
% matrix
[n_lips, ~, ~] = QuadratureData3D(Options.quadrature);% Local N of IPs
n_gips = n_lips * size(RVEMODEL.Conectivity,1); % Global N of IPs
int_vars = zeros(n_gips, 2);

end