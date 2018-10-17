%%% CO: facility location
clear all; close all; clc
ALG = 'CO_facloc';
addpath(genpath('..'));
rng(456);
m = 500;
n = 300;
C = sprandn(n, m, 0.01); % each column is a client location c_i
data.C = C;
data.m = m;
data.n = n;
x0 = randn(n*m, 1);
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'co-facloc');
param.mem_size = 5;
param.itermax = 500;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-8;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG, 6);





