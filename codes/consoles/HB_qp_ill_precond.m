%%% HB: random quadratic program / linear system
clear all; close all; clc
ALG = 'HB_qp_ill_precond';
addpath(genpath('..'));
rng(456);
n = 1000;
tol = 0.005;
B = randn(n/2, n);
A = B'*B + tol * eye(n, n);
b = randn(n, 1);
% diagonal scaling
D = diag(sum(abs(A),2));
A = D \ A;
b = D \ b;
E = diag(sum(abs(A),1));
A = A * E^(-1);
data.A = A;
data.b = b;
mu = tol;
L = norm(A, 'fro');
data.alpha = 4 / (sqrt(mu) + sqrt(L))^2;
data.beta = (sqrt(L) - sqrt(mu)) / (sqrt(L) + sqrt(mu));
data.n = n;
x0 = randn(2*n, 1);
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'hb-qp');
param.mem_size = 5;
param.itermax = 1000;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-5;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG, 2);





