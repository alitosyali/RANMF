clear
close all
clc

%% choose the data and corresponding parameters

[parentdir,~,~] = fileparts(pwd);

% % benchmark
% community = load(fullfile(parentdir, '/benchmark/community4.dat'));
% edges = load(fullfile(parentdir, '/benchmark/network4.dat'));
% true_label = community(:,2);
% N = size(community,1);
% A = zeros(N,N);
% for i = 1:size(edges,1)
%     A(edges(i,1),edges(i,2)) = 1;
% end
% r = length(unique(true_label));

% % pcn
% A = load(fullfile(parentdir, '/data/pcn_adj_mat.txt'));
% r = 20;

% % cornell
% A = load(fullfile(parentdir, '/data/cornell_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/cornell_labels.txt'));
% r = 5;

% % wisconsin
% A = load(fullfile(parentdir, '/data/wisconcin_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/wisconsin_labels.txt'));
% r = 5;

% % washington
% A = load(fullfile(parentdir, '/data/washington_adj_mat.txt'));
% true_label = load(fullfile(parentdir, '/data/washington_labels.txt'));
% r = 5;

% texas
A = load(fullfile(parentdir, '/data/texas_adj_mat.txt'));
true_label = load(fullfile(parentdir, '/data/texas_labels.txt'));
r = 5;

%% Algorithm

Kmax = r;
VV = GCSpectralClust1(A,Kmax);
T = [];
for cols = 2:size(VV,2)
    num_clusters = max(unique(VV(:,cols)));
    % for pcn
    %     db_score = db_index(A, VV(:,cols))
    %     quality_score = QFDistBased(VV(:,cols), A)
    
    % for others
    if max(unique(VV(:,cols))) == max(unique(true_label))
        predicted = bestMap(true_label, VV(:,cols));
        accuracy = sum(predicted == true_label)/length(predicted);
        Jaccard = PSJaccard(VV(:,cols), true_label);
        NMI = PSNMI(VV(:,cols), true_label);
        table(max(unique(VV(:,cols))),  Jaccard, NMI, accuracy)
    end
end

