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
plot(1:N, abs((xp - xp1)))
legend("abs(xp - xp1)")

fprintf("x: %d\t", max(log10(abs(xp - xp1))))
fprintf("G: %d\t", max(abs(G - G)))
fprintf("V: %d\t", max(abs(V - V1)))
fprintf("P: %d\t", max(abs(P - P1)))
fprintf("th: %d\n", max(abs(th - 1./th1)))



function [x, G, V, P, th, x1, G1, V1, P1, th1]=rkalman(A,B,C,D,y,c,tau)
%RKALMAN  robust Kalman estimator.
%
%   [x_e,G,V] = RKALMAN(A,B,C,D,y,c,tau) compute the "delayed" robust
%   Kalman estimator for the nominal discrete-time model
%
%      x[n+1] = Ax[n] + Bv[n]         {State equation}
%        y[n] = Cx[n] + Dv[n]         {Measurements}
%
%   with disturbance and measurement noise v with variance I. The robust
%   "delayed" estimator uses only past measurements up to y[n-1] to
%   generate the optimal robust estimate x_e[n] of x[n] and is easier to
%   embed in digital control loops. The equation of the robust "delayed"
%   estimator:
%
%      x_e[n+1|n] = Ax_e[n|n-1] + G (y[n] - Cx_e[n|n-1])
%
%   RKALMAN returns the delayed estimate x_e, the estimator gain G, the
%   least feavorable covariance matrix V of the estimation error
%   x[n+1]-x_e[n+1|n].
%
%   The robust estimator is designed knowing that the true model belongs to
%   a ball about the nominal one. The models in that ball are such that
%   the Tau-divergence between them and the nominal model is less than a
%   certain tolerance. To design the robust Kalman estimator, it is
%   necessary to specify:
%
%    - the tolerance c (striclty positive)
%    - the parameter tau (in the interval [0,1]) of the Tau-divergence
%    - the measurements y
%
%   For more details see: "Robust Kalman filtering under incremental
%   model perturbations" by M. Zorzi
%
%   See also RKITERATION, MAXTOL.
%   Author(s): Mattia Zorzi 20-8-2015
% check the inputs
if nargin==6
    tau=1;
end
% parameters
n=size(A,1);
p=size(C,1);
m=size(B,2);
T=size(y,1);
% 
sys.A = A;
sys.B = B;
sys.C = C;
sys.D = D;
sys.K = zeros(1,2);
% transform the model
% A=A-B*D'*(D*D')^-1*C;
% B=[(B*(eye(m)-D'*(D*D')^-1*D)*B')^0.5 zeros(n,m-n)];
% D=[zeros(p,n) (D*D')^0.5];
% Q=B*B';
% R=D*D';
% % check conditions
% assert(c>0,'Tolerance c must be positive')
% assert(rank(ctrb(A,Q))>=n,'The model must be reachable')
% assert(rank(ctrb(A',C'))>=n,'The model must be observable')
cN = maxtol(A,B,C,D,tau,2*n);
if c>cN
    warning('Tolerance c is too large: the filter gain may not exist')
end
% init
x=zeros(T,n);       x1=zeros(T,n);
V=zeros(n,n,T+1);   V1=zeros(n,n,T+1);
P=zeros(n,n,T+1);   P1=zeros(n,n,T+1);
G=zeros(n,p,T+1);   G1=zeros(n,p,T+1);
th=zeros(T,1);      th1=zeros(T,1);
V(:,:,1)=eye(n);    V1(:,:,1)=eye(n);
% iterative part
for k=1:T
    [x(k+1,:), V(:,:,k+1), G(:,:,k+1), P(:,:,k+1), th(k)]=rkiteration(A,B,C,D,V(:,:,k),0,c,x(k,:),y(k,:)); 
    [x1(k+1,:), V1(:,:,k+1), G1(:,:,k+1), P1(:,:,k+1), th1(k)]=RobustKalmanFilter(sys,V1(:,:,k),x1(k,:),y(k,:),c);
end
% resize
x=x(2:T+1,:);       x1=x1(2:T+1,:);
P=P(:,:,2:T+1);     P1=P1(:,:,2:T+1);
V=V(:,:,2:T+1);     V1=V1(:,:,2:T+1);
G=G(:,:,2:T+1);     G1=G1(:,:,2:T+1);
end