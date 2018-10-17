%%% PGD: convex-concave matrix games
clear all; close all; clc
ALG = 'PGD_ccm';
addpath(genpath('..'));
rng(456);
mp = 500;
np = 1500;
P = randn(mp, np);
e = ones(np, 1);
A = [P', eye(np, np), e];
n = mp + np + 1;
data.Pt = P';
data.P = P;
data.PPt = P*P';
data.Pe = P * e;
data.e = e;
data.alpha = 1.8 / max(eig(A'*A));
data.z = zeros(np, 1);
data.np = np;
data.mp = mp;
x0 = randn(n, 1);
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'pgd-ccm');
param.mem_size = 5;
param.itermax = 1000;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-5;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG, 2);


