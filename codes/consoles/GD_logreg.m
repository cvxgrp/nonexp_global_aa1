%%% GD: regularized logistic regression
clear all; close all; clc
ALG = 'GD_logreg';
addpath(genpath('..'));
rng(456);
data.xi = importdata('data/madelon_train.data.txt');
data.yi = importdata('data/madelon_train.labels.txt');
data.lambda = 1e-2;
n = size(data.xi, 2);
data.m = length(data.yi);
L = norm(data.xi)^2/4/data.m + data.lambda;
data.alpha = 2 / (L + data.lambda);
x0 = randn(n, 1);
x0 = x0 / norm(x0) * 1e-3; %zeros(n, 1);
F = @(x)fx(x,data,'gd-logreg');
param.mem_size = 5; 
param.itermax = 1000;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-5;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG);

%% plot objective values
m = length(data.yi);
fx = @(x)(sum(log(1+exp(-data.xi*x.*data.yi)))/m ...
    +data.lambda/2 * norm(x)^2/m);
fvalueplots(fx, param, x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, ALG);




