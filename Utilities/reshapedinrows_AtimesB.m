function Cvects = reshapedinrows_AtimesB( Avects , Bvects , dimsA , dimsB )

%auxiliar indexes:
indBaux   = reshape(reshape(1:prod(dimsB),dimsB)',1,prod(dimsB)); %col indexes of the transposed Bs
c_indexes = [  repmat((1:dimsA(1))',dimsB(2),1) , repelem((1:dimsB(2))',dimsA(1),1) ]; %[i,j]
%c_indexesA = repmat((1:dimsA(1))',dimsB(2),1) ;
%c_indexesB = repelem((1:dimsB(2))',dimsA(1),1) ; %[i,j]

%reshaping prior to performing the product:
%note that common dim = dimsA(2) = dimsB(1)
reshAmat = reshape(       Avects      , size(Avects,1) , dimsA(1) , dimsA(2) );
reshBmat = reshape( Bvects(:,indBaux) , size(Bvects,1) , dimsB(2) , dimsB(1) );

%compute coefficients of the product:
Cvects = sum( reshAmat(:,c_indexes(:,1),:) .* reshBmat(:,c_indexes(:,2),:) , 3 ); %preferred over than dot(), as supports one mat being an scalar

%Cvects = sum( repmat(reshAmat,1,dimsB(2),1) .* repelem(reshBmat,1,dimsA(1),1) , 3 );
%reshA  = reshAmat(:,c_indexesA,:) ;
%reshB  = reshBmat(:,c_indexesB,:);
%Cvects = sum(reshA.*reshB,3);
%Cvects = dot( reshAmat(:,c_indexes(:,1),:) , reshBmat(:,c_indexes(:,2),:) , 3 );

return