function [d, centroid] = compute_d(i, A, clusters)

% i: cluster index
% A: adjacency matrix
% clusters: set of clusters

this_cluster = clusters{i};
centroid = mean(A(this_cluster,:));
num_nodes_in_cluster = length(this_cluster);

d = 0;
for u = 1:num_nodes_in_cluster
    d = d + norm(A(this_cluster(u),:) - centroid);   
end

d = d/num_nodes_in_cluster;

end