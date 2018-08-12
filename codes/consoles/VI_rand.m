%% Value Iteration: random instances
clear all; close all; clc;
ALG = 'VI_rand';
addpath(genpath('/Users/junziz/Desktop/Programming/nonexp_aa/codes'));
rng(456);
data = struct();
S = 300; 
A = 200;
gamma = 0.99;
P = zeros(S, S, A);
for a = 1 : A
    tmp = sprand(S, S, 0.01);
    tmp = tmp + eye(S) * 1e-3; % ensure no NaN
    P(:, :, a) = diag(sum(tmp, 2))^(-1) * tmp;
end
R = sprandn(S, A, 0.01);
data.P = P;
data.R = R;
data.gamma = gamma;
x0 = randn(S, 1);
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'vi-rand');
param.mem_size = 5;
param.itermax = 50; 
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-5;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG);

