function thetay = phi(eta, theta)
%%% supplementary function for AA-I(-safe) updates (Powell regularization)
if abs(eta) >= theta
    thetay = 1;
elseif eta ~= 0
    thetay = (1-sign(eta)*theta) / (1-eta);
else
    thetay = 1 - theta;
end