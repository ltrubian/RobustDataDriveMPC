function [V_next, P_next, G, lambda] = riccati_iteration(sys, V, c)
% RICCATI_ITERATION Computes a single iteration of the robust Riccati equation.
%
%   [V_next, P_next, G, lambda] = riccati_iteration(sys, V, c)
%   Calculates the next least-favorable conditional covariance matrix and 
%   filter gain given the current covariance and the ambiguity set radius.
%
% INPUTS:
%   sys - State-space system structure containing A, B, C, D.
%   V   - Current least-favorable covariance matrix.
%   c   - Radius of the ambiguity set.
%
% OUTPUTS:
%   V_next - Next least-favorable conditional covariance matrix.
%   P_next - Next nominal conditional covariance matrix.
%   G      - Optimal filter gain matrix.
%   lambda - Lagrange multiplier associated with V and c, computed via 
%            a modified Regula-Falsi root-finding method.

arguments
    sys (1,1) struct
    V   (:,:) double
    c   (1,1) double {mustBeNonnegative(c)}
end

G = (sys.A*V*sys.C' + sys.B*sys.D') / (sys.C*V*sys.C' + sys.D*sys.D');

% nominal conditional covariance matrix
P_next = (sys.A-G*sys.C) * V * (sys.A-G*sys.C)' + (sys.B-G*sys.D) * (sys.B-G*sys.D)';

% START: COMPUTATION OF LAMBDA
% Lambda is the value such that
%   trace(inv(eye(n) - P/lambda) - eye(n)) + log(det(eye(n) - P/lambda)) = c*2
% This is computed using a variantion of the regula-falsi (similar to the
% secant method)
e = eig(P_next);
% the computation of the previous function can be simplified once P is
% diagonalized
gamfun = @(lam) sum(e./(lam-e) + log(1-e./lam)) - c*2;
% initialize values
lambda0 = max(e) * 1.01;    gamma0 = gamfun(lambda0);
lambda1 = lambda0 * 10;     gamma1 = gamfun(lambda1);
lambda = lambda0;

while abs(lambda0 - lambda1) >= 1e-9 && abs(gamma1 - gamma0) >= 1e-9
    lambda = (lambda0 * gamma1 - gamma0 * lambda1 ) / (gamma1 - gamma0);
    gamma2 = gamfun(lambda);
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
    lambda1 = lambda; gamma1 = gamma2;
end
% END: COMPUTATION OF LAMBDA

% least-favorable covariance matrix
V_next = (eye(size(sys.A)) - P_next/lambda) \ P_next;

% enforce symmetry and reduce numerical errors
V_next = (V_next+V_next')/2;
end