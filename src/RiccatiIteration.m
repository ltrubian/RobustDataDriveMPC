function [V_next, P_next, G, lambda] = RiccatiIteration(sys, V, c)
%RobustKalmanFilter one iteration
%
%   [V_next, G, P, lambda] = RobustKalmanFilter(sys, V, c)
%   compute one interation of the Kalman Filter
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

arguments
    sys (1,1) struct
    V   (:,:) double
    c   (1,1) double {mustBeNonnegative(c)}
end

G = (sys.A*V*sys.C' + sys.B*sys.D') / (sys.C*V*sys.C' + sys.D*sys.D');

% nominal conditional covariance matrix
P_next = (sys.A-G*sys.C)*V*(sys.A-G*sys.C)' + (sys.B-G*sys.D)*(sys.B-G*sys.D)';

% least-favorable covariance matrix
lambda = LagrangeMultiplier(P_next, c, "fast");
V_next = ( eye(size(sys.A)) - P_next/lambda ) \ P_next;

end