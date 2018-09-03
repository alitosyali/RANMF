function [S] = cosine_similarity(A)
    n = size(A,1); % number of nodes in given graph
    degree_matrix = sum(A,2);
    S = zeros(n,n);
    for node_i = 1:n
        for node_j = 1:n
            if (degree_matrix(node_i) * degree_matrix(node_j) ~= 0)
                S(node_i, node_j) = dot(A(node_i,:),A(node_j,:))/(norm(A(node_i,:))*norm(A(node_j,:)));
            end
        end
    end
end

