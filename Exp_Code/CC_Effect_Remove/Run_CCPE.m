addpath(genpath('drtoolbox'));
load('data/hNSC_preprocessed.mat')
lambda=70;
gamma=140;
sigma=0.001;  %Gaussian distribution
expression_matrix = double(expression_matrix)
[pseudotime] = CCPE(expression_matrix,lambda,gamma,sigma);
csvwrite("data/hNSC_pseudotime.csv", pseudotime);
