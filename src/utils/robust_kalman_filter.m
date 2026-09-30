function [x_pred, V_next, G, P_next, lambda] = robust_kalman_filter(sys, V, x, y, c)
% ROBUST_KALMAN_FILTER one iteration
%
%   [x_pred, V_next, G, P, lambda] = RobustKalmanFilter(sys, V, x_hat, y, c)
%   compute one prediction ahead of the Kalman Filter
%
% INPUT
%   sys:    struct with fields A, B, C, D, K
%   V:      least-favorable covariance matrix
%   x:      current state
%   y:      current output
%   c:      radius of the ambiguity set
%
% OUTPUT
%   x_pred: one-step ahead state prediction
%   V_next: next least-favorable conditional covariance matrix
%   G:      kalman gain used for prediction
%   P:      next nominal conditional covariance matrix
%   lambda: current lagrangian multiplier

arguments
    sys (1,1) struct
    V   (:,:) double
    x   (:,1) double
    y   (:,1) double
    c   (1,1) double {mustBeNonnegative(c)}
end

[V_next, P_next, G, lambda] = riccati_iteration(sys, V, c);

x_pred = sys.A*x + G * (y - sys.C*x);

end