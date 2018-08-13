function fvalueplots(fx, param, x_rec_origin, x_rec_aa1, x_rec_aa1_safe,...
    t_rec_origin, t_rec_aa1, t_rec_aa1_safe, ALG)
% plot objective values of the iterations for optimization problems
fx_rec_origin = zeros(param.itermax+1, 1);
fx_rec_aa1 = zeros(param.itermax+1, 1);
fx_rec_aa1_safe = zeros(param.itermax+1, 1);
for i = 1 : param.itermax+1
    fx_rec_origin(i) = fx(x_rec_origin(:,i));
    fx_rec_aa1(i) = fx(x_rec_aa1(:,i));
    fx_rec_aa1_safe(i) = fx(x_rec_aa1_safe(:,i));
end

% plots against iter
figure;
fstar = min([fx_rec_origin; fx_rec_aa1; fx_rec_aa1_safe]);
id_origin = min([param.itermax + 1, ...
    find([fx_rec_origin-fstar; 1e-11] < 1e-10, 1), ...
    find(isnan(fx_rec_origin) == 1, 1), ...
    find(abs(fx_rec_origin) == Inf, 1)]);
id_aa1 = min([param.itermax + 1, ...
    find([fx_rec_aa1-fstar; 1e-11] < 1e-10, 1), ...
    find(isnan(fx_rec_aa1) == 1, 1), ...
    find(abs(fx_rec_aa1) == Inf, 1)]);
id_aa1_safe = min([param.itermax + 1, ...
    find([fx_rec_aa1_safe-fstar; 1e-11] < 1e-10, 1), ...
    find(isnan(fx_rec_aa1_safe) == 1, 1), ...
    find(abs(fx_rec_aa1_safe) == Inf, 1)]);
semilogy(fx_rec_aa1(1:id_aa1) - fstar, 'b', 'LineWidth', 2); hold on
semilogy(fx_rec_aa1_safe(1:id_aa1_safe) - fstar, 'color', ...
[0.9100, 0.4100, 0.1700], 'LineWidth', 2);
semilogy(fx_rec_origin(1:id_origin) - fstar, 'r', 'LineWidth', 2);
xlim([0, param.itermax+1]);
ylim([10^(-8), 10^0]);
xlabel('iteration number', 'FontSize', 18);
ylabel('$F(x^k)-F(x^*)$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca,'YTick',10.^(-8:2:0));
legend('aa1', 'aa1-safe', 'original');
title('obj v.s. iter', 'FontSize', 18);
set(gca,'fontsize',18)
print('-dpdf', ['../figures/', ALG, '_iter_obj.pdf']);
hold off

% plots against time
figure;
semilogy(t_rec_aa1(1:id_aa1), ...
    fx_rec_aa1(1:id_aa1) - fstar, 'b', 'LineWidth', 2); hold on
semilogy(t_rec_aa1_safe(1:id_aa1_safe), ...
    fx_rec_aa1_safe(1:id_aa1_safe) - fstar, 'color', ...
[0.9100, 0.4100, 0.1700], 'LineWidth', 2);
semilogy(t_rec_origin(1:id_origin), ...
    fx_rec_origin(1:id_origin) - fstar, 'r', 'LineWidth', 2);
xlim([0, max([t_rec_origin(id_origin); t_rec_aa1(id_aa1); ...
    t_rec_aa1_safe(id_aa1_safe)])]);
ylim([10^(-8), 10^0]);
xlabel('time (seconds)', 'FontSize', 18);
ylabel('$F(x^k)-F(x^*)$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca,'YTick',10.^(-8:2:0));
legend('aa1', 'aa1-safe', 'original');
title('obj v.s. time', 'FontSize', 18);
set(gca,'fontsize',18)
print('-dpdf', ['../figures/', ALG, '_time_obj.pdf']);
hold off