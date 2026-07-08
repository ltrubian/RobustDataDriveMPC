function [x_pred, V_next, G, P_next, lambda] = RobustKalmanFilter(sys, V, x, y, u, c)
%RobustKalmanFilter one iteration
% 
%   [x_pred, V_next, G, P, lambda] = RobustKalmanFilter(sys, V, x_hat, y, c)
%   compute one prediction ahead of the Kalman Filter
%
% INPUT
%   sys:    struct with fields A, B, C, D, K
%   V:      least-favorable covariance matrix
%   x:      current state
%   y:      current output
%   u:      current input of the system
%   c:      radius of the ambiguity set
% OUTPUT
%   x_pred: prediction
%   V_next: next least-favorable conditional covariance matrix
%   P:      next nominal conditional covariance matrix
%   lambda: current lagrangian multiplier
% 
% DEV-STATUS: la funzione è stata testata: per 300 step mantiene la
% distanza dalla versione del professore con un errore di 1e-5
arguments
    sys (1,1) struct
    V   (:,:) double
    x   (:,1) double
    y   (:,1) double
    u   (:,1) double
    c   (1,1) double {mustBeNonnegative(c)}
end

[V_next, P_next, G, lambda] = RiccatiIteration(sys, V, c);

x_pred = sys.A*x + G * (y - sys.C*x) + sys.K * u;

end