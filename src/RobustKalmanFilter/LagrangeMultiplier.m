function lambda2 = LagrangeMultiplier(P, c, type)
%LagrangeMultiplier
%   lambda2 = LagrangeMultiplier(P, c, type) finds the lagrange multiplier
%   using secant method
%
%   type:   "normal" standard gamma function
%           "fast"   simplified version obtained after diagonalization of P
%                    (NOTE: that the eigenvalues are already computed to
%                    decide the starting point)
%
% DEV-STATUS: la funzione è stata testata nel senso che "fast" e "normal"
% riportano lo stesso risultato ed è lo zero della funzione (gamma - 2*c)
arguments
    P    (:,:) double
    c    (1,1) double {mustBeNonnegative(c)}
    type (1,:) string % {mustBeMember(type,["fast", "normal"])} = "normal"
end
n = size(P,1);
e = eig(P);

switch type
    case "normal"
        gamfun = @(lam) trace(inv(eye(n) - P/lam) - eye(n)) + log(det(eye(n) - P/lam)) - c*2;
    case "fast"
        gamfun = @(lam) sum(e./(lam-e)) + log(prod(1-e./lam)) - c*2;
end

lambda0 = max(e) * 1.01;
lambda1 = lambda0 * 10;
side = 0;
gamma0 = gamfun(lambda0);
gamma1 = gamfun(lambda1);
gamma2 = 1;
lambda2 = lambda0;
k = 0;

while abs(gamma2) >= 1e-9 && abs(lambda0 - lambda1) > 1e-9
    k = k + 1;
    lambda2 = (lambda0 * gamma1 - gamma0 * lambda1 ) / (gamma1 - gamma0);
    gamma2 = gamfun(lambda2);

    % --- next iteration values ---
    if sign(gamma2) == sign(gamma1)
        if side == -1
            mp = 1 - gamma2/gamma1;
            if mp < 0
                gamma0 = gamma0 * 0.5;
            else
                gamma0 = gamma0 * mp;
            end
        end
        lambda1 = lambda2;
        gamma1 = gamma2;
        side = -1;
    else
        if side == 1
            mp = 1 - gamma2/gamma0;
            if mp < 0
                gamma1 = gamma1 * 0.5;
            else
                gamma1 = gamma1 * mp;
            end
        end
        lambda0 = lambda2;
        gamma0 = gamma2;
        side = 1;
    end

end
end