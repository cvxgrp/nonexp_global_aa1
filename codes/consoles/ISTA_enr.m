%%% ISTA: elastic net regression
clear all; close all; clc
ALG = 'ISTA_enr';
addpath(genpath('..'));
rng(456);
q = 500;
p = 1000; 
A = randn(q, p);
xhat = sprandn(p, 1, 0.1);
b = A * xhat + randn(q, 1) * 0.1;
data.AtA = A' * A;
data.Atb = A' * b;
data.z = zeros(p, 1);
mumax = norm(A' * b, 'inf');
data.mu = 1e-3 * mumax;
L = max(eig(A'*A)) + data.mu/2;
data.alpha = 1.8 / L;
x0 = randn(p, 1);
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'ista-enr');
param.mem_size = 5;
param.itermax = 1000;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-8;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG, 2);






