function [ Stresses , Ctg ] = OrthotropicElastic3D_V( Strains , MaterialProps )

% MaterialProps Struct arrangement (all scalars):
% MaterialProps.Young1
% MaterialProps.Young2
% MaterialProps.Young3
% MaterialProps.Poisson23
% MaterialProps.Poisson12
% MaterialProps.Poisson13
% MaterialProps.ShearModuli13
% MaterialProps.ShearModuli12
% MaterialProps.ShearModuli23
% MaterialProps.AngleDegree


% Compute Intrinsec Ctg      
Ctg = [         -(- MaterialProps.Young3*MaterialProps.Young1^2*MaterialProps.Poisson23^2 + MaterialProps.Young2*MaterialProps.Young1^2)/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2), -(MaterialProps.Young2*(MaterialProps.Young1*MaterialProps.Young2*MaterialProps.Poisson12 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson13*MaterialProps.Poisson23))/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2), -(MaterialProps.Young2*MaterialProps.Young3*(MaterialProps.Young1*MaterialProps.Poisson13 + MaterialProps.Young1*MaterialProps.Poisson12*MaterialProps.Poisson23))/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2),   0,   0,   0;
  -(MaterialProps.Young1*MaterialProps.Poisson12*MaterialProps.Young2^2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson13*MaterialProps.Poisson23*MaterialProps.Young2)/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2),       -(MaterialProps.Young2*(- MaterialProps.Young2*MaterialProps.Young3*MaterialProps.Poisson13^2 + MaterialProps.Young1*MaterialProps.Young2))/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2), -(MaterialProps.Young2*MaterialProps.Young3*(MaterialProps.Young1*MaterialProps.Poisson23 + MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13))/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2),   0,   0,   0;
 -(MaterialProps.Young1*MaterialProps.Young2*MaterialProps.Young3*MaterialProps.Poisson13 + MaterialProps.Young1*MaterialProps.Young2*MaterialProps.Young3*MaterialProps.Poisson12*MaterialProps.Poisson23)/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2), -(MaterialProps.Young2*(MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23 + MaterialProps.Young2*MaterialProps.Young3*MaterialProps.Poisson12*MaterialProps.Poisson13))/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2),       -(MaterialProps.Young2*MaterialProps.Young3*(- MaterialProps.Young2*MaterialProps.Poisson12^2 + MaterialProps.Young1))/(MaterialProps.Young2^2*MaterialProps.Poisson12^2 + 2*MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson12*MaterialProps.Poisson13*MaterialProps.Poisson23 + MaterialProps.Young3*MaterialProps.Young2*MaterialProps.Poisson13^2 - MaterialProps.Young1*MaterialProps.Young2 + MaterialProps.Young1*MaterialProps.Young3*MaterialProps.Poisson23^2),   0,   0,   0;
                                                                                                                  0,                                                                                                                 0,                                                                                                              0, MaterialProps.ShearModuli12,   0,   0;
                                                                                                                  0,                                                                                                                 0,                                                                                                              0,   0, MaterialProps.ShearModuli13,   0;
                                                                                                                  0,                                                                                                                 0,                                                                                                              0,   0,   0, MaterialProps.ShearModuli23 ];

% Rotate Ctg
[ T , Tinv ] = ZRotationMatrices3D( MaterialProps.AngleDegree );
Ctg = Tinv * Ctg * diag([ 1 1 1 2 2 2 ]) * T * diag([ 1 1 1 0.5 0.5 0.5 ]);

% Repeat Ctgs along 3rd dimension
Ctg = reshape(Ctg,1,36) .* ones(size(Strains,1),1);

% Compute linear elastic stresses
Stresses = reshapedinrows_AtimesB( Ctg , Strains , [6 6] , [6 1]);

return
 