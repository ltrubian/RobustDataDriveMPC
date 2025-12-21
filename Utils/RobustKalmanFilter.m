function [x_pred, V_next, G, P, lambda] = RobustKalmanFilter(sys, V, x_hat, y, c)
% 
% 
% 
arguments
    sys     struct
    V       double
    x_hat   double
    y       double
    c       double
end

tmp = sys.C*V*sys.C' + sys.D*sys.D';

G = (sys.A*V*sys.C + sys.B*sys.D')' * inv(tmp);

x_pred = sys.A*x_hat + G * (y - sys.C*x_hat);

P = sys.A*V*sys.A' - G * tmp * G' + sys.B*sys.B';

lambda = LagrangeMultiplier(P, c);

V_next = inv(P - eye(size(A,1))/lambda );

end