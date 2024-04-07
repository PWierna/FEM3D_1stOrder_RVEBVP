posnodal3D = [-1 -1 -1; 1 -1 -1; 1  1 -1;-1  1 -1;...
                -1 -1  1; 1 -1  1; 1  1  1;-1  1  1;...
                 0 -1 -1; 1  0 -1; 0  1 -1;-1  0 -1;...
                -1 -1  0; 1 -1  0; 1  1  0;-1  1  0;...
                 0 -1  1; 1  0  1; 0  1  1;-1  0  1;...
                 0  0 -1; 0 -1  0; 1  0  0; 0  1  0;...
                -1  0  0; 0  0  1; 0  0  0 ];

face_normals_3D = [ -1,0,0 ; 1,0,0 ; 0,-1,0 ; 0,1,0 ; 0,0,-1 ; 0,0,1 ];

face_nodes_3D = zeros(6,9); %6 faces , 9 nodes per face
for f = 1:size(face_nodes_3D,1)
    face_normal = face_normals_3D(f,:);
    face_nodes = sort(find( posnodal3D(:,face_normal~=0) == face_normal(face_normal~=0) ));
    
    %Sort:
    face_nodes_coords = posnodal3D(face_nodes,face_normal==0);
    face_nodes_3D(f,1) = face_nodes( face_nodes_coords(
    
end


plot3(posnodal3D(:,1),posnodal3D(:,2),posnodal3D(:,3),'ok','LineWidth',1.5);hold on;grid on
xlabel('1');ylabel('2');zlabel('3');axis equal

for n=1:27
    text(posnodal3D(n,1),posnodal3D(n,2),posnodal3D(n,3),['\leftarrow ',char(string(n))],'FontSize',15)
end

f=1; %face
plot3(posnodal3D(face_nodes_3D(f,:),1),posnodal3D(face_nodes_3D(f,:),2),posnodal3D(face_nodes_3D(f,:),3),'or','LineWidth',1.5);
for fn=1:9
    text(posnodal3D(face_nodes_3D(f,fn),1),posnodal3D(face_nodes_3D(f,fn),2),posnodal3D(face_nodes_3D(f,fn),3),['\leftarrow ',char(string(face_nodes_3D(f,fn)))],'FontSize',15,'r')
end
