function y = fx(x, data, method)
%%% vanilla algorithms
if strcmp(method, 'gd-logreg')
    % gradient descent for regularized logistic regression
    xi = data.xi;
    yi = data.yi;
    alpha = data.alpha;
    lambda = data.lambda;
    m = data.m;
    tmp = xi * x .* yi;
    g = lambda * x + ((-yi ./ (exp(tmp) + 1))' * xi)'/m;
    y = x - alpha * g;
elseif strcmp(method, 'pgd-rand')
    % projected gradient descent for random non-negative least squares
    AA = data.AA;
    Ab = data.Ab;
    z = data.z;
    alpha = data.alpha;
    y = x - alpha * (AA * x - Ab);
    y = max([y, z], [], 2);
elseif strcmp(method, 'pgd-ccm')
    %projected gradient descent for convex-concave mat game
    Pt = data.Pt;
    P = data.P;
    PPt = data.PPt;
    Pe = data.Pe;
    np = data.np;
    mp = data.mp;
    e = data.e;
    z = data.z;
    alpha = data.alpha;
    u = x(1: mp);
    s = x(mp+1: mp+np);
    t = x(mp + np + 1);
    gu = PPt * u + P * s - Pe * t;
    gs = Pt * u + s - t * e;
    gt = np * t - Pe' * u - sum(s);
    u1 = u - alpha * gu;
    s1 = s - alpha * gs;
    t1 = t - alpha * gt;
    y = [SimplexProj(u1')'; max([s1, z], [], 2); t1];
elseif strcmp(method, 'ap-lp-sdhe')
    % alternating projection for LP as self-dual homogeneous embedding
    P = data.P;
    c = data.c;
    z = data.z;
    y = max([P * x + c, z], [], 2);
elseif strcmp(method, 'dap-lp-sdhe')
    % Dykstra alternating projection for LP as self-dual homogeneous embedding
    P = data.P;
    c = data.c;
    z = data.z;
    x12 = P * x + c;
    x1 = 2 * x12 - x;
    x1 = max([x1, z], [], 2);
    y = x + x1 - x12;
elseif strcmp(method, 'ista-enr')
    % ISTA for elastic net regression
    mu = data.mu;
    alpha = data.alpha;
    AtA = data.AtA;
    Atb = data.Atb;
    z = data.z;
    x1 = x - alpha * AtA * x + alpha * Atb - alpha*mu/2*x;
    y = sign(x1) .* max([abs(x1) - alpha*mu/2, z], [], 2);
elseif strcmp(method, 'admm-huber')
    % ADMM for robust (Huber) regression
    X = data.X;
    yi = data.y;
    Rup = data.R_up;
    Q = data.Q;
    n = data.n;
    m = data.m;
    e = data.e;
    theta = x(1:n);
    u = x(n+1:n+m);
    % update of t
    t1tilde = yi - X * theta - u;
    t1 = t1tilde - t1tilde .* min([e/2,1./abs(t1tilde)], [], 2);
    % update of theta
    Qyut = Q' * (yi - u - t1);
    if m >= n
        theta1 = Rup \ Qyut(1:n);
    else
        theta1 = Rup \ Qyut;
    end
    % update of u
    u1 = u + t1 + X * theta1 - yi;
    y = [theta1; u1];
elseif strcmp(method,'scs-lp-sdhe')
    % SCS for LP as self-dual homoegeneous embedding
    IQinv = data.IQinv;
    n = data.n;
    z = data.z;
    u0 = x(1:n);
    v0 = x(n+1:2*n);
    utilde1 = IQinv * (u0+v0);
    u1 = max([utilde1 - v0, z], [], 2);
    v1 = v0 - utilde1 + u1;
    y = [u1; v1];
elseif strcmp(method,'scs-socp-sdhe')
    % SCS for SOCP as self-dual homoegeneous embedding
    IQinv = data.IQinv;
    n = data.n;
    n0 = data.n0;
    m0 = data.m0;
    u0 = x(1:n);
    v0 = x(n+1:2*n);
    utilde1 = IQinv * (u0+v0);
    utilde1_v0 = utilde1 - v0;
    u1 = [utilde1_v0(1:n0); soc_proj(utilde1_v0(n0+1:n0+m0)); ...
        max([utilde1_v0(n0+m0+1), 0])];
    v1 = v0 - utilde1 + u1;
    y = [u1; v1];
elseif strcmp(method, 'vi-rand')
    % value iteration for random MDP
    P = data.P;
    R = data.R;
    gamma = data.gamma;
    y = max(R + gamma * ...
        squeeze(sum(bsxfun(@times, P, reshape(x, 1, length(x))), 2)),[],2);
elseif strcmp(method, 'hb-qp')
    % heavy-ball for simple QP / linear systems
    A = data.A;
    b = data.b;
    alpha = data.alpha;
    beta = data.beta;
    n = data.n;
    x0 = x(1:n);
    x1_0 = x(n+1:end);
    x1 = x1_0 - alpha * (A * x1_0 + b) + beta * (x1_0 - x0);
    x0 = x1_0;
    y = [x0; x1];
elseif strcmp(method, 'co-facloc')
    % consensus optimization for facility location problem
    C = data.C;
    m = data.m;
    n = data.n;
    X = reshape(x, [n, m]); % reshape so that each column is z_i
    X1 = bsxfun(@prox_norm, X + C, zeros(1,m)) - C;
    X2 = X + 2*repmat(mean(X1,2), [1,m]) - X1 - repmat(mean(X,2), [1,m]);
    y = reshape(X2, [n*m, 1]);
end