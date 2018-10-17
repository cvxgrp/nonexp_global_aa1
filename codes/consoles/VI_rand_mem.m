%% Value Iteration: random instances -- memory tests
clear all; close all; clc;
ALG = 'VI_rand_mem';
addpath(genpath('..'));
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
param.itermax = 50;
res0 = norm(x0 - F(x0));
mems = [2, 5, 10, 20, 50];

%% algorithm comparisons
tol = 1e-5;
% Algorithms & Residual computation
param_origin = param;
param_origin.itermax = param.itermax * 2;
[x_rec_origin, t_rec_origin] = alg_iter(x0, F, param_origin, 'original');
res_origin = zeros(param_origin.itermax+1, 1);
res_aa1 = cell(length(mems), 1);
res_aa1_safe = cell(length(mems), 1);
t_rec_aa1s = cell(length(mems), 1);
t_rec_aa1_safes = cell(length(mems), 1);

% original algorithm
count_origin = 1;
for i = 1 : param_origin.itermax+1
    res_origin(i) = norm(x_rec_origin(:,i) - F(x_rec_origin(:,i)));
    count_origin = count_origin + 1;
    if res_origin(i) < tol * res0 ...
            || isnan(res_origin(i)) || res_origin(i)==Inf
        count_origin = count_origin - 1;
        break
    end
end
t_rec_origin0 = t_rec_origin;
res_origin = res_origin(1:min(count_origin, param_origin.itermax+1));
t_rec_origin = t_rec_origin(1:length(res_origin));

% AA-I & AA-I-safe
param.theta = 0.01;
param.tau = 0.001;
param.D = 1e6;
param.epsilon = 1e-6;
for mem = 1 : length(mems)
    res_aa1{mem} = zeros(param.itermax+1, 1);
    res_aa1_safe{mem} = zeros(param.itermax+1, 1);
    param.mem_size = mems(mem);
    [x_rec_aa1, t_rec_aa1] = alg_iter(x0, F, param, 'aa1');
    [x_rec_aa1_safe, t_rec_aa1_safe, rec_aa1_safe] ...
        = alg_iter(x0, F, param, 'aa1-safe');
    count_aa1 = 1;
    for i = 1 : param.itermax+1
        res_aa1{mem}(i) = norm(x_rec_aa1(:,i) - F(x_rec_aa1(:,i)));
        count_aa1 = count_aa1 + 1;
        if res_aa1{mem}(i) < tol * res0 ...
                || isnan(res_aa1{mem}(i)) || res_aa1{mem}(i)==Inf
            count_aa1 = count_aa1 - 1;
            break
        end
    end
    count_aa1_safe = 1;
    for i = 1 : param.itermax+1
        res_aa1_safe{mem}(i) = ...
            norm(x_rec_aa1_safe(:,i) - F(x_rec_aa1_safe(:,i)));
        count_aa1_safe = count_aa1_safe + 1;
        if res_aa1_safe{mem}(i) < tol * res0 ...
                || isnan(res_aa1_safe{mem}(i)) || res_aa1_safe{mem}(i)==Inf
            count_aa1_safe = count_aa1_safe - 1;
            break
        end
    end
    t_rec_aa10 = t_rec_aa1;
    t_rec_aa1_safe0 = t_rec_aa1_safe;
    res_aa1{mem} = res_aa1{mem}(1:min(count_aa1, param.itermax+1));
    res_aa1_safe{mem} = ...
        res_aa1_safe{mem}(1:min(count_aa1_safe, param.itermax+1));
    t_rec_aa1s{mem} = t_rec_aa1(1:length(res_aa1{mem}));
    t_rec_aa1_safes{mem} = t_rec_aa1_safe(1:length(res_aa1_safe{mem}));
end

% plot the original curve only to max(itermax, time-aa1, time-aa1-safe)
t_max = 0;
for mem = 1 : length(mems)
    t_max = max([t_max; t_rec_aa1s{mem}; t_rec_aa1_safes{mem}]);
end
count_origin = max(sum(t_rec_origin <= t_max), param.itermax+1);
res_origin = res_origin(1:count_origin);
t_rec_origin = t_rec_origin(1:count_origin);


%% Plots
% res against iter
colors = ['b', 'm', 'k', 'y', 'g'];
figure;
semilogy(res_origin(1:param.itermax+1)/res0, 'r', 'LineWidth', 2); hold on
for i = 1 : length(mems)
    if i ~= 2
        semilogy(res_aa1{i}/res0, [colors(i), '-.'], 'LineWidth', 2)
        semilogy(res_aa1_safe{i}/res0, [colors(i), '-'], 'LineWidth', 2);
    else
        semilogy(res_aa1{i}/res0, 'color', [0.9100, 0.4100, 0.1700], ...
            'LineStyle', '-.', 'LineWidth', 2)
        semilogy(res_aa1_safe{i}/res0, 'color', [0.9100, 0.4100, 0.1700], ...
            'LineStyle', '-', 'LineWidth', 2);
    end
end
xlabel('iteration number', 'FontSize', 18);
ylabel('$\|g(x^k)\|_2/\|g(x^0)\|_2$', 'Interpreter', 'latex', ...
    'FontSize', 18);
legend({'original', 'aa1,m=2', 'aa1-safe, m=2', ...
    'aa1,m=5', 'aa1-safe, m=5', ...
    'aa1,m=10', 'aa1-safe, m=10', ...
    'aa1,m=20', 'aa1-safe, m=20', ...
    'aa1,m=50', 'aa1-safe, m=50'}, 'FontSize', 12);
title(['res v.s. iter, res0=', num2str(res0, '%3.2e')], 'FontSize', 18);
set(gca, 'fontsize', 18);
print('-dpdf', ['../figures/', ALG, '_iter.pdf']);
hold off

% res against time
figure;
semilogy(t_rec_origin, res_origin/res0, 'r', 'LineWidth', 2); hold on
xlimup = max(t_rec_origin);
for i = 1 : length(mems)
    if i ~= 2
        semilogy(t_rec_aa1s{i}, res_aa1{i}/res0, [colors(i), '-.'], ...
            'LineWidth', 2);
        semilogy(t_rec_aa1_safes{i}, res_aa1_safe{i}/res0, ...
            [colors(i), '-'], 'LineWidth', 2);
        xlimup = max([xlimup; t_rec_aa1s{i}; t_rec_aa1_safes{i}]);
    else
        semilogy(t_rec_aa1s{i}, res_aa1{i}/res0, 'color', ...
            [0.9100, 0.4100, 0.1700], 'LineStyle', '-.', 'LineWidth', 2);
        semilogy(t_rec_aa1_safes{i}, res_aa1_safe{i}/res0, 'color', ...
            [0.9100, 0.4100, 0.1700], 'LineStyle', '-', 'LineWidth', 2);
        xlimup = max([xlimup; t_rec_aa1s{i}; t_rec_aa1_safes{i}]);
    end
end
xlim([0, xlimup]);
xlabel('time (seconds)', 'FontSize', 18);
ylabel('$\|g(x^k)\|_2/\|g(x^0)\|_2$', 'Interpreter', 'latex', ...
    'FontSize', 18);
legend({'original', 'aa1,m=2', 'aa1-safe, m=2', ...
    'aa1,m=5', 'aa1-safe, m=5', ...
    'aa1,m=10', 'aa1-safe, m=10', ...
    'aa1,m=20', 'aa1-safe, m=20', ...
    'aa1,m=50', 'aa1-safe, m=50'}, 'FontSize', 12);
title(['res v.s. time, res0=', num2str(res0, '%3.2e')], 'FontSize', 18);
set(gca,'fontsize',18)
print('-dpdf', ['../figures/', ALG, '_time.pdf']);
hold off

