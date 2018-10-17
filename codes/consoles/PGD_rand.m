%%% PGD: random non-negative least squares
clear all; close all; clc
ALG = 'PGD_rand';
addpath(genpath('..'));
rng(456);
m = 500;
n = 1000;
r = 10;
A = randn(m, n);
b = randn(m, 1);
data.AA = A' * A;
data.Ab = A' * b;
data.alpha = 1.8 / max(eig(data.AA));
data.z = zeros(n, 1);
x0 = randn(n, 1);
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'pgd-rand');
param.mem_size = 5;
param.itermax = 1000;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-5;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG, 2);





