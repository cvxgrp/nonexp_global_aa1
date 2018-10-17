function [x_rec_origin, x_rec_aa1, x_rec_aa1_safe, ...
    t_rec_origin0, t_rec_aa10, t_rec_aa1_safe0, rec_aa1_safe] ...
    = console_common(tol, x0, F, param, res0, ALG, ratio_origin)
%% Algorithms
% ratio_origin (>=1): run the original algorithm for ratio_origin times 
% the iteration numbers of the aa1 and aa1-safe algorithms (useful only for
% plotting purposes, otherwise just set to 1)
param_origin = param;
param_origin.itermax = param.itermax * ratio_origin;
[x_rec_origin, t_rec_origin] = alg_iter(x0, F, param_origin, 'original');
[x_rec_aa1, t_rec_aa1] = alg_iter(x0, F, param, 'aa1');
param.theta = 0.01;
param.tau = 0.001;
param.D = 1e6;
param.epsilon = 1e-6;
[x_rec_aa1_safe, t_rec_aa1_safe, rec_aa1_safe] ...
    = alg_iter(x0, F, param, 'aa1-safe');
res_origin = zeros(param_origin.itermax+1, 1);
res_aa1 = zeros(param.itermax+1, 1);
res_aa1_safe = zeros(param.itermax+1, 1);

%% Residual computation
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
count_aa1 = 1;
for i = 1 : param.itermax+1
    res_aa1(i) = norm(x_rec_aa1(:,i) - F(x_rec_aa1(:,i)));
    count_aa1 = count_aa1 + 1;
    if res_aa1(i) < tol * res0 ...
            || isnan(res_aa1(i)) || res_aa1(i)==Inf
        count_aa1 = count_aa1 - 1;
        break
    end
end
count_aa1_safe = 1;
for i = 1 : param.itermax+1
    res_aa1_safe(i) = norm(x_rec_aa1_safe(:,i) - F(x_rec_aa1_safe(:,i)));
    count_aa1_safe = count_aa1_safe + 1;
    if res_aa1_safe(i) < tol * res0 ...
            || isnan(res_aa1_safe(i)) || res_aa1_safe(i)==Inf
        count_aa1_safe = count_aa1_safe - 1;
        break
    end
end

t_rec_origin0 = t_rec_origin(1:param.itermax+1);
t_rec_aa10 = t_rec_aa1;
t_rec_aa1_safe0 = t_rec_aa1_safe;

res_origin = res_origin(1:min(count_origin, param_origin.itermax+1));
res_aa1 = res_aa1(1:min(count_aa1, param.itermax+1));
res_aa1_safe = res_aa1_safe(1:min(count_aa1_safe, param.itermax+1));
t_rec_origin = t_rec_origin(1:length(res_origin));
t_rec_aa1 = t_rec_aa1(1:length(res_aa1));
t_rec_aa1_safe = t_rec_aa1_safe(1:length(res_aa1_safe));

% plot the original curve only to max(itermax, time-aa1, time-aa1-safe)
t_max = max(t_rec_aa1(end), t_rec_aa1_safe(end));
count_origin = max(sum(t_rec_origin <= t_max), param.itermax+1);
res_origin = res_origin(1:count_origin);
t_rec_origin = t_rec_origin(1:count_origin);


%% Plots
% res against iter
figure;
semilogy(res_aa1/res0, 'b', 'LineWidth', 2); hold on
semilogy(res_aa1_safe/res0, 'color', ...
[0.9100, 0.4100, 0.1700], 'LineWidth', 2);
semilogy(res_origin(1:param.itermax+1)/res0, 'r', 'LineWidth', 2);
xlim([0, param.itermax+1]);
ylim_min = min([res_aa1; res_aa1_safe; res_origin] / res0)/2;
ylim_max = max([res_aa1; res_aa1_safe; res_origin] / res0)*2;
xlabel('iteration number', 'FontSize', 18);
ylabel('$\|g(x^k)\|_2/\|g(x^0)\|_2$', 'Interpreter', 'latex', ...
     'FontSize', 18);
exp_max = min(ceil(log(ylim_max) / log(10)), 5);
exp_min = max(floor(log(ylim_min) / log(10)), -16);
ylim([10^exp_min, 10^exp_max])
set(gca,'YTick',10.^(exp_min:2:exp_max));
legend('aa1', 'aa1-safe', 'original');
title(['res v.s. iter, res0=', num2str(res0, '%3.2e')], 'FontSize', 18);
set(gca,'fontsize',18)
print('-dpdf', ['../figures/', ALG, '_iter.pdf']);
hold off

% res against time
time_ratio_aa1 = t_rec_aa10(end) / t_rec_origin0(end);
time_ratio_aa1_safe = t_rec_aa1_safe0(end) / t_rec_origin0(end);
figure;
semilogy(t_rec_aa1, res_aa1/res0, 'b', 'LineWidth', 2); hold on
semilogy(t_rec_aa1_safe, res_aa1_safe/res0, 'color', ...
[0.9100, 0.4100, 0.1700], 'LineWidth', 2);
semilogy(t_rec_origin, res_origin/res0, 'r', 'LineWidth', 2);
xlim([0, max([t_rec_origin; t_rec_aa1; t_rec_aa1_safe])]);
ylim_min = min([res_aa1; res_aa1_safe; res_origin] / res0)/2;
ylim_max = max([res_aa1; res_aa1_safe; res_origin] / res0)*2;
ylim([ylim_min, ylim_max])
xlabel('time (seconds)', 'FontSize', 18);
ylabel('$\|g(x^k)\|_2/\|g(x^0)\|_2$', 'Interpreter', 'latex', ...
    'FontSize', 18);
exp_max = min(ceil(log(ylim_max) / log(10)), 5);
exp_min = max(floor(log(ylim_min) / log(10)), -16);
ylim([10^exp_min, 10^exp_max])
set(gca,'YTick',10.^(exp_min:2:exp_max));
legend('aa1', 'aa1-safe', 'original');
title({['res v.s. time, res0=', num2str(res0, '%3.2e')], ...
    ['time ratio: aa1 = ', num2str(time_ratio_aa1, '%3.2e'), ...
    ', aa1-safe = ', num2str(time_ratio_aa1_safe, '%3.2e')]}, ...
    'FontSize', 18);
set(gca,'fontsize',18)
print('-dpdf', ['../figures/', ALG, '_time.pdf']);
hold off