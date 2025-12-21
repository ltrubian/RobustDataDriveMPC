function lambda = LagrangeMultiplier(P, c)
n = size(P,1);
l_prec = 1;
lambda = 2;
gamma_prec = trace(inv(eye(n) - P/l_prec) - eye(n)) + log(det(eye(n)) )
while abs() >= 1e-9

end


end