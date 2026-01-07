function nodes12sorted = pairup_nodes( nodes1 , nodes2 , coordinates12 )

% Retrieve coordinates for both node sets:
coords1 = coordinates12(nodes1,:);
coords2 = coordinates12(nodes2,:);

% Initialize sorted indexes for nodes2:
nodes2sorted = zeros(length(nodes2),1);

% Sort by rows:
for i = 1:length(nodes2)
    
    % Retrieve coordinates of master node:
    masterc_i = coords1(i,:);
    
    % Find distance w.r.t. all the slave nodes:
    dist2master = vecnorm(coords2-masterc_i,2,2);

    % Identify point with the minimum distance:
    slave_i = find( dist2master == min(dist2master) );
    
    % Retrieve:
    if isscalar(slave_i==1)
        nodes2sorted(i) = nodes2( slave_i ); 
    else
        error('More than 1 match found for a single master node');
    end
    
end

% Retrieve sorted node pairs:
nodes12sorted = [nodes1 nodes2sorted];



end