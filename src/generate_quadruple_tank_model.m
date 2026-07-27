% generate_quadruple_tank_model.m
% Generates the discrete-time state-space matrices A, B, C, D, K 
% for the Quadruple-Tank Process (Johansson, 2000).
% 
% The resulting model matches the structure:
% x(t+1) = A x(t) + B v(t) + K u(t)
% y(t)   = C x(t) + D v(t)
% where the joint noise covariance [B; D] is square and invertible.

clear; clc;

%% Configuration
is_minimum_phase = false; % Set to false for non-minimum phase
Ts = 1.0; % Sampling time in seconds

%% Physical Parameters (Johansson, 2000)
A1 = 28; A3 = 28; % Area of tanks 1 and 3 [cm^2]
A2 = 32; A4 = 32; % Area of tanks 2 and 4 [cm^2]
a1 = 0.071; a3 = 0.071; % Cross-section of outlet hole 1 and 3 [cm^2]
a2 = 0.057; a4 = 0.057; % Cross-section of outlet hole 2 and 4 [cm^2]
kc = 0.50; % Sensor calibration [V/cm]
g = 981; % Gravitational constant [cm/s^2]
k1 = 3.33; % Pump 1 constant [cm^3/Vs]
k2 = 3.35; % Pump 2 constant [cm^3/Vs]

%% Operating Points
if is_minimum_phase
    disp('Generating Minimum Phase Quadruple-Tank Model...');
    gamma1 = 0.70;
    gamma2 = 0.60;
    h10 = 12.4; h20 = 12.7; h30 = 1.8; h40 = 1.4;
else
    disp('Generating Non-Minimum Phase Quadruple-Tank Model...');
    gamma1 = 0.43;
    gamma2 = 0.34;
    h10 = 12.6; h20 = 13.0; h30 = 4.8; h40 = 4.9;
end

%% Continuous-Time Linearized Model
% Time constants
T1 = (A1/a1) * sqrt(2*h10/g);
T2 = (A2/a2) * sqrt(2*h20/g);
T3 = (A3/a3) * sqrt(2*h30/g);
T4 = (A4/a4) * sqrt(2*h40/g);

Ac = [-1/T1,      0, A3/(A1*T3),          0;
          0, -1/T2,          0, A4/(A2*T4);
          0,      0,     -1/T3,          0;
          0,      0,          0,     -1/T4];

Bc = [gamma1*k1/A1,            0;
                 0, gamma2*k2/A2;
                 0, (1-gamma2)*k2/A3;
      (1-gamma1)*k1/A4,            0];

Cc = [kc, 0, 0, 0;
       0, kc, 0, 0];

Dc = zeros(2, 2);

%% Discretization
% Convert to discrete time using zero-order hold
sys_c = ss(Ac, Bc, Cc, Dc);
sys_d = c2d(sys_c, Ts, 'zoh');

A = sys_d.A;
K = sys_d.B; % K is the input-to-state matrix
C = sys_d.C;

%% Noise Matrices (B and D)
% We synthesize a full-rank joint covariance matrix for process and measurement noise.
% This guarantees that [B; D] is square and invertible.
n = 4; % States
p = 2; % Outputs

% Define arbitrary small variances for the physical process and sensors
% (e.g. 1e-4 for process noise on levels, 1e-3 for sensor noise)
Q = diag([1e-4, 1e-4, 1e-4, 1e-4]); 
R = diag([1e-3, 1e-3]);

% Joint covariance matrix
Sigma = blkdiag(Q, R);

% Cholesky factorization to find B and D
NoiseMat = chol(Sigma, 'lower');

B = NoiseMat(1:n, :);
D = NoiseMat(n+1:end, :);

%% Display and Save
disp('---------------------------------------------------------');
disp(['Matrix [B; D] size: ', num2str(size(NoiseMat, 1)), ' x ', num2str(size(NoiseMat, 2))]);
disp(['Rank of [B; D]: ', num2str(rank(NoiseMat))]);
disp('---------------------------------------------------------');

if is_minimum_phase
    save('quadruple_tank_min_phase.mat', 'A', 'K', 'B', 'C', 'D', 'Sigma', 'Ts');
    disp('Model saved to quadruple_tank_min_phase.mat');
else
    save('quadruple_tank_nonmin_phase.mat', 'A', 'K', 'B', 'C', 'D', 'Sigma', 'Ts');
    disp('Model saved to quadruple_tank_nonmin_phase.mat');
end
