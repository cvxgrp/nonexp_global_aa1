function y = prox_norm(x, emp)
% proximal operator of the l2 norm
% emp is not used: only for the sake of bsxfun
y = max(1-1/norm(x), 0) * x + emp - emp;
end