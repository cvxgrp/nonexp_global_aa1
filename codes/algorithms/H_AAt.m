function y = H_AAt(x, Hv1, Hv2)
%%% supplementary function for AA-I(-safe) updates
if isempty(Hv1) || isempty(Hv2)
    y = x;
else
    y = x + Hv2 * (Hv1' * x);
end
end