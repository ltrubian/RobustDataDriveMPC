function [V_next, P_next, G, lambda] = RiccatiIteration(sys, V, c)
%RobustKalmanFilter one iteration
%
%   [V_next, P_next, G, lambda] = RiccatiIteration(sys, V, c)
%   compute one interation of the Riccati equation
%
% INPUT
%   sys:    struct with fields A, B, C, D
%   V:      least-favorable covariance matrix
%   c:      radius of the ambiguity set
%
% OUTPUT
%   V_next: next least-favorable conditional covariance matrix
%   G:      filter gain
%   P_next: next nominal conditional covariance matrix
%   lambda: lagrange multiplier associated with V and c

arguments
    sys (1,1) struct
    V   (:,:) double
    c   (1,1) double {mustBeNonnegative(c)}
end

G = (sys.A*V*sys.C' + sys.B*sys.D') / (sys.C*V*sys.C' + sys.D*sys.D');

% nominal conditional covariance matrix
P_next = (sys.A-G*sys.C)*V*(sys.A-G*sys.C)' + (sys.B-G*sys.D)*(sys.B-G*sys.D)';

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
V_next = ( eye(size(sys.A)) - P_next/lambda ) \ P_next;

% enforce symmetry and reduce numerical errors
V_next = (V_next+V_next')/2;
end