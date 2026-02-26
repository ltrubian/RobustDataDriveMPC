function [x_pred, V_next, G, P, lambda] = RobustKalmanFilter(sys, V, x, y, c)
%RobustKalmanFilter one iteration
% 
%   [x_pred, V_next, G, P, lambda] = RobustKalmanFilter(sys, V, x_hat, y, c)
%   compute one interation of the Kalman Filter
%
% INPUT
%   sys:    struct with fields A, B, C, D
%   V:      least-favorable covariance matrix
%   x:      current state
%   y:      current output
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
    x   (1,:) double
    y   (1,:) double
    c   (1,1) double {mustBePositive(c)}
end

tmp = sys.C*V*sys.C' + sys.D*sys.D';

G = (sys.A*V*sys.C + sys.B*sys.D')' / (tmp);

x_pred = sys.A*x + G * (y - sys.C*x);

% nominal conditional covariance matrix
P = sys.A*V*sys.A' - G * tmp * G' + sys.B*sys.B';

% least-favorable covariance matrix
lambda = LagrangeMultiplier(P, c, "fast");
V_next = inv( inv(P)-eye(size(sys.A))/lambda );

end