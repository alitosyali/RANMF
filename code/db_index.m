function [DB] = db_index(A, predict_label)

cluster_labels = unique(predict_label);
num_clusters = length(cluster_labels);

clusters = [];
for k = 1:num_clusters
    clusters = [clusters {find(predict_label==cluster_labels(k))}];
end

D_ij = [];

for i = 1:num_clusters
    [d_i, centroid_i] = compute_d(i, A, clusters);
    for j = i+1:num_clusters
        [d_j, centroid_j] = compute_d(j, A, clusters);
        if norm(centroid_i - centroid_j) ~= 0            
            D_ij = [D_ij (d_i + d_j)/norm(centroid_i - centroid_j)];
        end
    end
end

DB = (1/num_clusters) * max(D_ij);

end




