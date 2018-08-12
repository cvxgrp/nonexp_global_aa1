function y = soc_proj(x)
% projection to second order cones
n = length(x);
x1 = x(1:n-1);
xn = x(n);
if norm(x1) <= xn
    y = x;
elseif norm(x1) <= -xn
    y = zeros(n, 1);
else
    y = (norm(x1) + xn) / 2 * [x1/norm(x1); 1];
end