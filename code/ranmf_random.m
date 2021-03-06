clear
close all
clc

%% choose the data and corresponding parameters

[parentdir,~,~] = fileparts(pwd);

% % benchmark
% community = load(fullfile(parentdir, 'benchmark/community3.dat'));
% edges = load(fullfile(parentdir, 'benchmark/network3.dat'));
% true_label = community(:,2);
% N = size(community,1);
% A = zeros(N,N);
% for i = 1:size(edges,1)
%     A(edges(i,1),edges(i,2)) = 1;
% end
% r = length(unique(true_label));
% num_iter = 20;
% num_runs = 20;
% lambda = 0.1;

% % pcn
% A = load(fullfile(parentdir, '/data/pcn_adj_mat.txt'));
% r = 15;
% lambda = 0.1;
% num_iter = 100;
% num_runs = 100;

% % cornell
% A = load(fullfile(parentdir, '/data/cornell_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/cornell_labels.txt'));
% num_iter = 100;
% num_runs = 100;
% lambda = 10;
% r = 5;

% % wisconsin
% A = load(fullfile(parentdir, '/data/wisconsin_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/wisconsin_labels.txt'));
% num_iter = 100;
% num_runs = 100;
% lambda = 25*50;
% r = 5;

% % washington
% A = load(fullfile(parentdir, '/data/washington_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/washington_labels.txt'));
% num_iter = 100;
% num_runs = 100;
% lambda = 1;
% r = 5;

% % texas
% A = load(fullfile(parentdir, '/data/texas_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/texas_labels.txt'));
% num_iter = 100;
% num_runs = 100;
% lambda = 3220;
% r = 5;

tic;
% number of nodes
n = size(A,1);

%% similarity matrix

% cosine
[S] = cosine_similarity(A);

% % katz
% beta = 0.1;
% [S] = citation_similarity(A, beta);

% % adjacency
% S = A + A' + eye(n);

%% updating the objective function

% Graph Laplacian
D = diag(sum(S,2));
L = D - S;

objective = zeros(num_runs,num_iter);
quality_score = zeros(1,num_runs);
% quality_score = zeros(3,num_runs);

for run = 1:num_runs
    fprintf('----------- Run: %d ----------\n', run)
    X_initial = abs(rand(n,r));
    S_initial = abs(rand(r,r));
    for iter = 1:num_iter
        if mod(iter,2) == 1
            numer_s = (X_initial')*(A)*(X_initial);
            denom_s = (X_initial')*(X_initial)*(S_initial)*(X_initial')*(X_initial);
            for k = 1:r
                for j = 1:r
                    S_initial(k,j) = S_initial(k,j)*(numer_s(k,j)/max(denom_s(k,j),realmin));
                end
            end
        else
            numer_x = (A)*(X_initial)*(S_initial') + (A')*(X_initial)*(S_initial) + lambda*(S')*X_initial;
            denom_x = (X_initial)*(S_initial)*(X_initial')*(X_initial)*(S_initial') + (X_initial)*(S_initial')*(X_initial')*(X_initial)*(S_initial) + 2*lambda*(D')*X_initial;
            for i = 1:n
                for k = 1:r
                    X_initial(i,k) = X_initial(i,k)*(numer_x(i,k)/max(denom_x(i,k),realmin))^(0.25);
                end
            end
        end
        objective(run,iter) = norm(A-(X_initial)*(S_initial)*(X_initial)')^2 + lambda*trace((X_initial')*(L)*(X_initial));
    end
    
    E = diag(sum(X_initial,1));
    W = X_initial/E;
    [~,predict_label] = max(W,[],2);
    
%     quality_score(run) = QFDistBased(predict_label, A);
%     quality_score(run) = db_index(A, predict_label);
    
        quality_score(1, run) = PSJaccard(predict_label, true_label); % jaccard
        if length(unique(predict_label)) == length(unique(true_label))
            quality_score(2, run) = PSNMI(predict_label, true_label); % nmi
        end
        predicted = bestMap(true_label, predict_label);
        quality_score(3, run) = sum(predicted == true_label)/length(predicted); % accuracy
end
toc;

% % for pcn
% mean(quality_score)
% table(mean(quality_score), 'RowNames', {'RANMF_Rand'})

% for others
jaccard = mean(quality_score(1, :));
nmi = mean(nonzeros(quality_score(2, :)));
accuracy = mean(quality_score(3, :));
table(jaccard, nmi, accuracy, 'RowNames', {'RANMF'})











