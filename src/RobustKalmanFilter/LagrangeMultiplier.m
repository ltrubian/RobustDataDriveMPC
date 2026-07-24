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
% DEV-STATUS: The function has been tested: "fast" and "normal" modes
% yield the same result, which is the root of the function (gamma - 2*c)
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
        gamfun = @(lam) sum(e./(lam-e) + log(1-e./lam)) - c*2;
end

lambda0 = max(e) * 1.01;
lambda1 = lambda0 * 10;
gamma0 = gamfun(lambda0);
gamma1 = gamfun(lambda1);
gamma2 = 1;
lambda2 = lambda0;

while abs(gamma2) >= 1e-9 && abs(lambda0 - lambda1) > 1e-9

    if abs(gamma1 - gamma0) < 1e-12
        break;
    end
    lambda2 = (lambda0 * gamma1 - gamma0 * lambda1 ) / (gamma1 - gamma0);
    gamma2 = gamfun(lambda2);

    % --- next iteration values ---
    if sign(gamma2) ~= sign(gamma1)
        lambda0 = lambda1; gamma0 = gamma1;
    else
        % regula-falsi, Anderson & Björk modification
        m = 1 - gamma2/gamma1;
        if m > 0
            gamma0 = gamma0 * m;
        else
            gamma0 = gamma0 / 2;
        end
    end
    lambda1 = lambda2; gamma1 = gamma2;
end
end