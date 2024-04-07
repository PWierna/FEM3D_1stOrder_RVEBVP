function StrainsM_gips  = compute_macro_strains_at_ips( GStrainM , GIPCoords , Orders )

%INPUTS:
%Orders    => [ order x , order y , order z ]
%GIPCoords => Rows=IP number ; Cols=[ X Y Z ]
%GStrainM  => [ Membrane(xx,yy,xy) , Bending(xx,yy,xy) , Shear(xz,yz) ]
%             [      (1,2,3)       ,      (4,5,6)      ,    (7,8)     ]

%OUTPUT:
%StrainsM_gips => Rows=IP number ; Cols=Strain[ xx , yy , zz , xy , xz , yz ]
StrainsM_gips      =    zeros( size(GStrainM,1) , 6 );
StrainsM_gips(:,1) =    GStrainM(1) - GIPCoords(:,3)*GStrainM(4);  %Exx
StrainsM_gips(:,2) =    GStrainM(2) - GIPCoords(:,3)*GStrainM(5);  %Eyy
StrainsM_gips(:,4) = 2*(GStrainM(3) - GIPCoords(:,3)*GStrainM(6)); %Exy 
%StrainsM_gips(:,3) = (zeros);                                     %Ezz
StrainsM_gips(:,5) =    GStrainM(7) - GIPCoords(:,1)*GStrainM(4) - GIPCoords(:,2)*GStrainM(6); %Exz
StrainsM_gips(:,6) =    GStrainM(8) - GIPCoords(:,1)*GStrainM(6) - GIPCoords(:,2)*GStrainM(5); %Eyz 



end



