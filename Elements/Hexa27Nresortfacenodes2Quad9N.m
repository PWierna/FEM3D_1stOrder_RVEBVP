function fidx_resort2D = Hexa27Nresortfacenodes2Quad9N

posnodal3D = [-1 -1 -1; 1 -1 -1; 1  1 -1;-1  1 -1;...
                -1 -1  1; 1 -1  1; 1  1  1;-1  1  1;...
                 0 -1 -1; 1  0 -1; 0  1 -1;-1  0 -1;...
                -1 -1  0; 1 -1  0; 1  1  0;-1  1  0;...
                 0 -1  1; 1  0  1; 0  1  1;-1  0  1;...
                 0  0 -1; 0 -1  0; 1  0  0; 0  1  0;...
                -1  0  0; 0  0  1; 0  0  0 ];

face_normals_3D = [ -1,0,0 ; 1,0,0 ; 0,-1,0 ; 0,1,0 ; 0,0,-1 ; 0,0,1 ];
%                     odd     even    even    odd      odd     even

face_connects_3D = zeros(6,9); %6 faces , 9 nodes per face
for f = 1:size(face_connects_3D,1)
    face_normal = face_normals_3D(f,:);
    face_connects_3D(f,:) = sort(find( posnodal3D(:,face_normal~=0) == face_normal(face_normal~=0) ));
end

% [1 4 5 8 12 13 16 20 25;
%                             2 3 6 7 10 14 15 18 23;
%                             1 2 5 6 9 13 14 17 22;
%                             3 4 7 8 11 15 16 19 24;
%                             1 2 3 4 9 10 11 12 21;
%                             5 6 7 8 17 18 19 20 26 ];
            
posnodal2D_even = [ -1 -1 ; 0 -1 ; 1 -1 ; 1  0 ; 1  1 ; 0  1 ; -1  1 ; -1  0 ; 0  0 ];
%posnodal2D_odd  = [ 1 -1 ; 0 -1 ; -1 -1 ; -1  0 ; -1  1 ; 0  1 ; 1  1 ; 1  0 ; 0  0 ];
posnodal2D_odd  = [ -1 -1 ; 0 -1 ; 1 -1 ; 1  0 ; 1  1 ; 0  1 ; -1  1 ; -1  0 ; 0  0 ];

evenfaces_idx = [ 2 , 3 , 6 ];
oddfaces_idx  = [ 1 , 4 , 5 ];

%compute resorting indices:
fidx_resort2D = zeros( length(face_normals_3D) , size(posnodal2D_even,1) );

% for f = 1:length(face_normals_3D) %for each 3D face
%     
%     fnormal  = face_normals_3D(f,:);
%     fconnect = face_connects_3D(f,:);
%     fcoords  = posnodal3D( fconnect , : );
%     fcoords  = fcoords( : , fnormal==0 ); %Get rid of the coord constant over the face
% 
%     for n = 1:size(posnodal2D_even,1) %for every node of this face
%         fidx_resort2D(f,n) = find(prod( fcoords == posnodal2D_even(n,:) , 2 ));
%     end
%     
% end
    
for f = evenfaces_idx %for each "even" 3D face
    
    fnormal  = face_normals_3D(f,:);
    fconnect = face_connects_3D(f,:);
    fcoords  = posnodal3D( fconnect , : );
    fcoords  = fcoords( : , fnormal==0 ); %Get rid of the coord constant over the face
    for n = 1:size(posnodal2D_even,1) %for every node of this face
        fidx_resort2D(f,n) = find(prod( fcoords == posnodal2D_even(n,:) , 2 ));
    end
    
end  

for f = oddfaces_idx %for each "odd" 3D face
    
    fnormal  = face_normals_3D(f,:);
    fconnect = face_connects_3D(f,:);
    fcoords  = posnodal3D( fconnect , : );
    fcoords  = fcoords( : , fnormal==0 ); %Get rid of the coord constant over the face
    for n = 1:size(posnodal2D_odd,1) %for every node of this face
        fidx_resort2D(f,n) = find(prod( fcoords == posnodal2D_odd(n,:) , 2 ));
    end
    
end

% resorted_face_connects_3D = face_connects_3D;
% for f=1:6
%     resorted_face_connects_3D(f,:) = face_connects_3D(f,fidx_resort2D(f,:)) -1 ;
% end

return