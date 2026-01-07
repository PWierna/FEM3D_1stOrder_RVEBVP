function nodepairsAB = pairupFaceNodes_bysorting( nodesA , nodesB , coordinates12 )

% Hard-setting: Round-off number of digits
rounddigits = 12;

%Get 1-2 coordinates for both faces:
faceA_ncoords = coordinates12( nodesA , : );
faceB_ncoords = coordinates12( nodesB , : );

%Force unique coordinate values by rounding-off floating point differencies:
faceA_ncoords = round(faceA_ncoords,rounddigits);
faceB_ncoords = round(faceB_ncoords,rounddigits);

%Sort by 1st-2nd coordinates:
[ ~ , idxsortA ] = sortrows( faceA_ncoords );
[ ~ , idxsortB ] = sortrows( faceB_ncoords );

%Retrieve paired nodes([ face , xmax face ]):
nodepairsAB = [ nodesA(idxsortA) , nodesB(idxsortB) ];

end