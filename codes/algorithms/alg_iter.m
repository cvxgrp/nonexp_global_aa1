function [x_rec, t_rec, rec] = alg_iter(x0, F, param, algorithm)
%%% acceleration algorithm iterations
% 'original': original algorithm
% 'aa1': AA-I
% 'aa1-safe': AA-I-safe
% 'aa2': AA-II
% 'aa2-reg': AA-II with regularization
beta = 0.1; % biased towards more progressive iterations
itermax = param.itermax;
n = length(x0);
x_rec = zeros(n, itermax+1);
x_rec(:, 1) = x0;
t_rec = zeros(itermax+1, 1);
rec = struct();

if strcmp(algorithm, 'original')
    t0 = cputime;
    t_rec(1) = 0;
    for i = 1 : itermax
        if mod(i, 20) == 1
            fprintf('iteration = %d\n', i);
        end
        x0 = F(x0);
        x_rec(:, i+1) = x0;
        t_rec(i+1) = cputime - t0;
    end
    
elseif strcmp(algorithm, 'aa1')
    mem_size = param.mem_size;
    Smem = [];
    Ymem = [];
    g0 = x0 - F(x0);
    t0 = cputime;
    t_rec(1) = 0;
    for i = 1 : itermax
        if mod(i, 20) == 1
            fprintf('iteration = %d\n', i);
        end
        if i == 1
            x1 = F(x0);
        else
            x1 = x0 - g0 - (Smem - Ymem) * ((Smem'*Ymem) \ (Smem' * g0));
        end
        if i <= mem_size - 1
            Smem(:, i) = x1 - x0;
            Ymem(:, i) = x1 - F(x1) - g0;
        else
            Smem = [Smem(:, 2:end), x1 - x0];
            Ymem = [Ymem(:, 2:end), x1 - F(x1) - g0];
        end
        x_rec(:, i+1) = x1;
        t_rec(i+1) = cputime - t0;
        x0 = x1;
        g0 = x0 - F(x0);
        
    end
    
elseif strcmp(algorithm, 'aa1-safe')
    mem_size = param.mem_size;
    theta = param.theta;
    tau = param.tau;
    D = param.D;
    epsilon = param.epsilon;
    g0 = x0 - F(x0);
    Ubar = norm(g0);
    Shat_mem = [];
    H_vecs1 = [];
    H_vecs2 = [];
    count = 0;
    m = 0;
    t0 = cputime;
    t_rec(1) = 0;
    rec.restart = [];
    rec.safeguard = [];
    for i = 1 : itermax
        m = m + 1; % default increase of memory
        if mod(i, 20) == 1
            fprintf('###iteration = %d\n', i);
        end
        if i == 1
            x1 = F(x0);
        else
            x1 = x0 - H_AA(g0, H_vecs1, H_vecs2);
        end
        
        %%% Safeguard checking
        Ubar0 = norm(g0);
        if Ubar0 <= D * Ubar * (count+1)^(-1-epsilon) || i == 1
            count = count + 1;
        else
            rec.safeguard = [rec.safeguard, i];
            x1 = beta * x0 + (1-beta) * F(x0);
        end
        x_rec(:, i+1) = x1;
        t_rec(i+1) = cputime - t0;
        s0 = x1 - x0;
        g1 = x1 - F(x1);
        y0 = g1 - g0;
        x0 = x1;
        g0 = x0 - F(x0);
        
        %%% Restart checking
        if m <= mem_size
            if ~isempty(Shat_mem) % only do when nonempty memory
                s0hat = s0 - ((Shat_mem * s0)' * Shat_mem)';
            else
                fprintf('iter=%d: no orthogonalization\n', i);
                s0hat = s0;
            end
            if norm(s0hat) < tau * norm(s0) % restart
                rec.restart = [rec.restart, i];
                s0hat = s0;
                m = 1;
                Shat_mem = [];
                H_vecs1 = [];
                H_vecs2 = [];
            end
        else % memory exceeds
            s0hat = s0;
            m = 1;
            Shat_mem = [];
            H_vecs1 = [];
            H_vecs2 = [];
        end
        Shat_mem = [Shat_mem; s0hat' / norm(s0hat)];
        
        %%% Powell regularization
        gamma0 = s0hat' * H_AA(y0, H_vecs1, H_vecs2) / norm(s0hat)^2;
        theta0 = phi(gamma0, theta); 
        y0tilde = theta0 * y0 - (1-theta0) * g0;
        
        %%% Update H_vecs
        Hytilde = H_AA(y0tilde, H_vecs1, H_vecs2);
        hvec1 = s0 - Hytilde;
        hvec2 = H_AAt(s0hat, H_vecs1, H_vecs2) / (s0hat' * Hytilde);
        H_vecs1 = [H_vecs1, hvec1];
        H_vecs2 = [H_vecs2, hvec2];
    end
    
elseif strcmp(algorithm, 'aa2')
    mem_size = param.mem_size;
    mem = x0;
    Fmem = F(x0);
    t0 = cputime;
    t_rec(1) = 0;
    for i = 1 : itermax
        if mod(i, 20) == 1
            fprintf('iteration = %d\n', i);
            param.print = true;
        else
            param.print = false;
        end
        Fmem_mem = Fmem - mem;
        e = ones(size(Fmem_mem,2), 1);
        alpha = (Fmem_mem' * Fmem_mem) \ e;
        alp = alpha / sum(alpha);
        x0 = beta * mem * alp + (1-beta) * Fmem * alp;
        if i <= mem_size - 1
            mem(:, i+1) = x0;
            Fmem(:, i+1) = F(x0);
        else
            mem = [mem(:, 2:end), x0];
            Fmem = [Fmem(:, 2:end), F(x0)];
        end
        x_rec(:, i+1) = x0;
        t_rec(i+1) = cputime - t0;
    end
    
elseif strcmp(algorithm, 'aa2-reg')
    mem_size = param.mem_size;
    mu = param.mu;
    mem = x0;
    Fmem = F(x0);
    t0 = cputime;
    t_rec(1) = 0;
    for i = 1 : itermax
        if mod(i, 20) == 1
            fprintf('iteration = %d\n', i);
            param.print = true;
        else
            param.print = false;
        end
        Fmem_mem = Fmem - mem;
        e = ones(size(Fmem_mem,2), 1);
        en = e; en(1) = 0; en = diag(en);
        alpha = (Fmem_mem' * Fmem_mem + mu * en) \ e;
        alp = alpha / sum(alpha);
        x0 = beta * mem * alp + (1-beta) * Fmem * alp;
        if i <= mem_size - 1
            mem(:, i+1) = x0;
            Fmem(:, i+1) = F(x0);
        else
            mem = [mem(:, 2:end), x0];
            Fmem = [Fmem(:, 2:end), F(x0)];
        end
        x_rec(:, i+1) = x0;
        t_rec(i+1) = cputime - t0;
    end
end