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
% DEV-STATUS: la funzione non è stata testata in nessun modo
arguments
    sys struct
    V   double
    x   double
    y   double
    c   double
end

tmp = sys.C*V*sys.C' + sys.D*sys.D';

G = (sys.A*V*sys.C + sys.B*sys.D')' / (tmp);

x_pred = sys.A*x + G * (y - sys.C*x);

% nominal conditional covariance matrix
P = sys.A*V*sys.A' - G * tmp * G' + sys.B*sys.B';

% least-favorable covariance matrix
lambda = LagrangeMultiplier(P, c);
V_next = inv(P - eye(size(A,1))/lambda );

end