%%% AP: LP by SDHE
clear all; close all; clc;
ALG = 'AP_lp_sdhe';
addpath(genpath('/Users/junziz/Desktop/Programming/nonexp_aa/codes'));
rng(456);
m0 = 500; n0 = 1000;
A0 = sprandn(m0, n0, 0.1);
zstar = randn(n0, 1);
z0 = zeros(n0, 1);
xstar = max([zstar, z0], [], 2);
sstar = max([-zstar, z0], [], 2);
ystar = randn(m0, 1);
b0 = A0 * xstar;
c0 = A0' * ystar + sstar;
% diagonal scaling of A, b and c
D = diag(sum(abs(A0), 2));
A0 = D \ A0;
b0 = D \ b0;
E = diag(sum(abs(A0), 1));
A0 = A0 * E^(-1);
c0 = E \ c0;
zn0 = zeros(n0, 1);
zm0 = zeros(m0, 1);
Zmm0 = zeros(m0, m0);
Zmn0 = zeros(m0, n0);
Znn0 = zeros(n0, n0);
en0 = ones(n0, 1);
I0 = eye(n0, n0);
Atilde = [Znn0, -A0', c0, -I0, zn0; ...
    A0, Zmm0, -b0, Zmn0, zm0; ...
    -c0', b0', 0, zn0', -1];
data = struct();
n = 2 * n0 + m0 + 2;
m = m0 + n0 + 2;
btilde = zeros(m-1, 1);
I = eye(n, n);
P = I - Atilde' * ((Atilde * Atilde') \ Atilde);
c = Atilde' * ((Atilde * Atilde') \ btilde);
z = zeros(n, 1);
z(n0+1:n0+m0) = -Inf;
data.P = P;
data.c = c;
data.z = z;
x0 = randn(n, 1); 
x0 = x0 / norm(x0);
F = @(x)fx(x,data,'ap-lp-sdhe');
param.mem_size = 5;
param.itermax = 1000;
res0 = norm(x0 - F(x0));

%% algorithm comparisons
tol = 1e-5;
[x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG);






