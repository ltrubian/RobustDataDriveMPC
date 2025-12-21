function gamma = gammafun(lam, P, c, type)
n = size(P,1);
e = eig(P);
if strcmp(type, "normal")
    gamma = trace(inv(eye(n) - P/lam) - eye(n)) + log(det(eye(n) - P/lam)) -2*c;
else
    gamma = sum(e./(lam-e)) + log(prod(1-e/lam)) - 2*c;
end
end