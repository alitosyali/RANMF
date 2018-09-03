clear
close all
clc

%% choose the data and corresponding parameters

[parentdir,~,~] = fileparts(pwd);

% % benchmark
% community = load(fullfile(parentdir, '\benchmark\community3.dat'));
% edges = load(fullfile(parentdir, '\benchmark\network3.dat'));
% true_label = community(:,2);
% N = size(community,1);
% A = zeros(N,N);
% for i = 1:size(edges,1)
%     A(edges(i,1),edges(i,2)) = 1;
% end
% r = length(unique(true_label));
% num_iter = 200;
% num_runs = 20;

% % pcn
% A = load(fullfile(parentdir, '/data/pcn_adj_mat.txt'));
% r = 20;
% num_iter = 100;
% num_runs = 100;

% % cornell
% A = load(fullfile(parentdir, '/data/cornell_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/cornell_labels.txt'));
% num_iter = 100;
% num_runs = 100;
% r = 5;

% % wisconsin
% A = load(fullfile(parentdir, '/data/wisconcin_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/wisconsin_labels.txt'));
% num_iter = 100;
% num_runs = 100;
% r = 5;

% % washington
% A = load(fullfile(parentdir, '/data/washington_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/washington_labels.txt'));
% num_iter = 100;
% num_runs = 100;
% r = 5;

% texas
A = load(fullfile(parentdir, '/data/texas_adj_mat.txt'));
true_label = load(fullfile(parentdir, '/data/texas_labels.txt'));
num_iter = 100;
num_runs = 100;
r = 5;

%% iterative algorithm

% number of nodes in given network
n = size(A,1);

objective = zeros(num_runs,num_iter);
% quality_score = zeros(1,num_runs);
quality_score = zeros(3,num_runs);
tic;
for run = 1:num_runs
    fprintf('----------- Run: %d ----------\n', run)
    X = abs(rand(n,r));
    S = abs(rand(r,r));
    for iter = 1:num_iter
        if mod(iter,2) == 1
            numer_s = (X')*(A)*(X);
            denom_s = (X')*(X)*(S)*(X')*(X);
            for k = 1:r
                for j = 1:r
                    S(k,j) = S(k,j)*(numer_s(k,j)/max(denom_s(k,j),realmin));
                end
            end
        else
            numer_x = (A)*(X)*(S') + (A')*(X)*(S);
            denom_x = (X)*(S)*(X')*(X)*(S') + (X)*(S')*(X')*(X)*(S);
            for i = 1:n
                for k = 1:r
                    X(i,k) = X(i,k)*(numer_x(i,k)/max(denom_x(i,k),realmin))^(0.25);
                end
            end
        end
        objective(run,iter) = norm(A-(X)*(S)*(X)')^2;
    end
    E = diag(sum(X,1));
    W = X/E;
    [~,predict_label] = max(W,[],2);
    
    if length(unique(predict_label)) == length(unique(true_label))
        %     quality_score(run) = QFDistBased(predict_label, A);
        %     quality_score(run) = db_index(A, predict_label);
        
        quality_score(1, run) = PSJaccard(predict_label, true_label); % jaccard
        quality_score(2, run) = PSNMI(predict_label, true_label); % nmi
        predicted = bestMap(true_label, predict_label);
        quality_score(3, run) = sum(predicted == true_label)/length(predicted); % accuracy
    end
end
toc;
% % for pcn
% db_score = mean(quality_score);
% table(db_score, 'RowNames', {'ANMF'})

% for others
jaccard = mean(quality_score(1, :));
nmi = mean(quality_score(2, :));
accuracy = mean(quality_score(3, :));
table(jaccard, nmi, accuracy, 'RowNames', {'ANMF'})
