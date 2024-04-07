function [ T , Tinv ] = ZRotationMatrices3D( Theta )

%Direction cosines:
%m  = (cos(+Theta*pi/180));
%n  = (sin(+Theta*pi/180));
%mi = (cos(-Theta*pi/180));
%ni = (sin(-Theta*pi/180));


%Compute transformation matrices (rotation about Z axis):
T = [      cos(Theta*pi/180)^2                     sin(Theta*pi/180)^2              0      2*cos(Theta*pi/180)*sin(Theta*pi/180)              0                       0             ; 
           sin(Theta*pi/180)^2                     cos(Theta*pi/180)^2              0     -2*cos(Theta*pi/180)*sin(Theta*pi/180)              0                       0             ;  
                   0                                        0                       1                         0                               0                       0             ;
      -cos(Theta*pi/180)*sin(Theta*pi/180)    cos(Theta*pi/180)*sin(Theta*pi/180)   0     cos(Theta*pi/180)^2-sin(Theta*pi/180)^2             0                       0             ;
                     0                                        0                     0                         0                          cos(Theta*pi/180)     sin(Theta*pi/180)    ;
                     0                                        0                     0                         0                         -sin(Theta*pi/180)     cos(Theta*pi/180)    ];


Tinv = [      cos(Theta*pi/180)^2                     sin(Theta*pi/180)^2              0     -2*cos(Theta*pi/180)*sin(Theta*pi/180)              0                       0             ; 
              sin(Theta*pi/180)^2                     cos(Theta*pi/180)^2              0      2*cos(Theta*pi/180)*sin(Theta*pi/180)              0                       0             ;  
                      0                                        0                       1                         0                               0                       0             ;
          cos(Theta*pi/180)*sin(Theta*pi/180)   -cos(Theta*pi/180)*sin(Theta*pi/180)   0     cos(Theta*pi/180)^2-sin(Theta*pi/180)^2             0                       0             ;
                      0                                        0                     0                         0                          cos(Theta*pi/180)    -sin(Theta*pi/180)    ;
                      0                                        0                     0                         0                          sin(Theta*pi/180)     cos(Theta*pi/180)    ];

end