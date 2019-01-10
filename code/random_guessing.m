clear;
close all;
clc;

% set the seed for reproducibility 
rng(42);

%% choose the dataset
[parentdir,~,~] = fileparts(pwd);

% benchmark
community = load(fullfile(parentdir, '/benchmark/community3.dat'));
edges = load(fullfile(parentdir, '/benchmark/network3.dat'));
true_label = community(:,2);
N = size(community,1);
A = zeros(N,N);
for i = 1:size(edges,1)
    A(edges(i,1),edges(i,2)) = 1;
end
r = length(unique(true_label));

% % pcn
% A = load(fullfile(parentdir, '/data/pcn_adj_mat.txt'));
% r = 15;

% % cornell
% A = load(fullfile(parentdir, '/data/cornell_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/cornell_labels.txt'));
% r = length(unique(true_label));

% % wisconsin
% A = load(fullfile(parentdir, '/data/wisconsin_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/wisconsin_labels.txt'));
% r = length(unique(true_label));

% % washington
% A = load(fullfile(parentdir, '/data/washington_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/washington_labels.txt'));
% r = length(unique(true_label));

% % texas
% A = load(fullfile(parentdir, '/data/texas_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/texas_labels.txt'));
% r = length(unique(true_label));

%% random guessing
N = size(A,1);
num_runs = 100;
% all_results = zeros(num_runs, 2); % for pcn
all_results = zeros(num_runs, 3); % for other datasets

for this_run = 1:num_runs
    
    predict_label = randi(r, N, 1);
    
%     % for PCN
%     all_results(this_run, :) = [QFDistBased(predict_label, A) db_index(A, predict_label)];
    
    % for others
    jaccard = PSJaccard(predict_label, true_label); % jaccard
    nmi = PSNMI(predict_label, true_label); % nmi
    predicted = bestMap(true_label, predict_label);
    accuracy = sum(predicted == true_label)/length(predicted); % accuracy
    all_results(this_run, :) = [jaccard, nmi, accuracy];
end

mean(all_results)

