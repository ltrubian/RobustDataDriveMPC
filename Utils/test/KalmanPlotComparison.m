%--------------------------------------------------------------------------
% Name:            RKF_Demo.m
%
% Description:     Application of the robust Kalman filter to the task of
%                  tracking a sampled Wiener process on noisy measurements.
%
% Author:          Mattia Zorzi
%
% Date:            Agoust 20, 2015
%--------------------------------------------------------------------------
addpath("../")
%%  True model
rng(30)
% Process noise variance
Q = 0.8;
% Measurement noise variance
R = 2;
%% Data generation
% Number of iterations
N = 300;
% True State and measurement
x = zeros(N,1);
y = zeros(N,1);
% True initial state
x(1) = randn;
% First measurement
y(1) = x(1) + sqrt(Q)*randn;
% Update true state and measurements
for i=2:N
    x(i) = x(i-1) + sqrt(Q)*randn;
    y(i) = x(i) + sqrt(R)*randn;
end
%% Nominal noise/process covariance
% Process noise variance
Qi = 2.0;
% Measurement noise variance
Ri = 1.6;
%% Robust Kalman filtering
% Initial apriori state estimate
% profile on
[xp, G, V, P, th, xp1, G1, V1, P1, th1] = rkalman(1,[sqrt(Qi) 0],1,[0 sqrt(Ri)],y,1e-10,0);
% profile viewer
%% Plot Results
figure
% Plot true state
b = plot(1:N,x(1:N));
hold on
% Plot estimates
r = plot(2:N,xp(1:N-1));
plot(2:N,xp1(1:N-1));
% Plot measurements
g = plot(1:N,y(1:N), '+');
title(['Robust Kalman Filtering: Tracking a sampled Wiener process']);
legend('True Value','RKF Estimate','RKF Estimate 2','Measurements');
xlabel('Time')
ylabel('Value')
grid on

figure(Name="lambda")
plot(1:N, th, 1:N, 1./th1)
legend("prof", "myfilter")

figure(Name="prediction difference")
plot(1:N, (abs(xp - xp1)))
legend("abs(xp - xp1)")

fprintf("x: %d\t", max(log10(abs(xp - xp1))))
fprintf("G: %d\t", max(abs(G - G)))
fprintf("V: %d\t", max(abs(V - V1)))
fprintf("P: %d\t", max(abs(P - P1)))
fprintf("th: %d\n", max(abs(th - 1./th1)))