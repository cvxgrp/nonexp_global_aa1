%%% SCS: SOCP by SDHE
clear all; close all; clc;
ALG = 'SCS_socp_sdhe';
addpath(genpath('..'));
rng(456);
m0 = 500; n0 = 700;
n = m0 + n0 + 1;
A = [sprandn(m0, n0/2, 0.1), eye(m0, n0/2)] + 1e-3 * randn(m0, n0);
zstar = randn(m0, 1);
sstar = soc_proj(zstar);
ystar = sstar - zstar;
xstar = randn(n0, 1);
b = A * xstar + sstar;
c = -A' * ystar;
% diagonal scaling of A, b and c
D = diag(sum(abs(A), 2));
A = D \ A;
b = D \ b;
E = diag(sum(abs(A), 1));
A = A * E^(-1);
c = E \ c;
Znn = zeros(n0, n0);
Zmm = zeros(m0, m0);
Q = [Znn, A', c; -A, Zmm, b; -c', -b', 0]; 
IQinv = inv(eye(n)+Q);
data.IQinv = IQinv;
data.n = n;
data.n0 = n0;
data.m0 = m0;
x0 = randn(2*n, 1);
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'scs-socp-sdhe');
param.mem_size = 5; 
param.itermax = 1000;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-5;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG, 2);

