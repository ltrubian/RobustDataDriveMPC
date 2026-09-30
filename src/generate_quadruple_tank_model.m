function generate_quadruple_tank_model()
% GENERATE_QUADRUPLE_TANK_MODEL Generates and saves the discrete-time Quadruple-Tank model.
%
% Constructs the state-space matrices for the Quadruple-Tank Process. 
% It linearizes the continuous-time dynamics around a defined operating 
% point and discretizes the system using a zero-order hold.
%
% Additionally, it synthesizes process (B) and measurement (D) noise gain 
% matrices using Cholesky factorization on a randomly generated positive 
% definite covariance matrix. This ensures the joint noise matrix [B; D] 
% is square and invertible.
%
% The resulting state-space matrices (A, K, B, C, D), covariance (Sigma), 
% and sample time (Ts) are saved to 'quadruple_tank.mat'.

%% Configuration
Ts = 1; % sampling time in seconds

% physical Parameters
A1 = 28; A3 = 28; % area of tanks 1 and 3 [cm^2]
A2 = 32; A4 = 32; % area of tanks 2 and 4 [cm^2]
a1 = 0.071; a3 = 0.071; % cross-section of outlet hole 1 and 3 [cm^2]
a2 = 0.057; a4 = 0.057; % cross-section of outlet hole 2 and 4 [cm^2]
g = 981; % gravitational constant [cm/s^2]

% operating Point
gamma1 = 0.70; gamma2 = 0.60;
h10 = 12.4; h20 = 12.7; h30 = 1.8; h40 = 1.4; % initial conditions [cm]
k1 = 3.33; % pump 1 constant [cm^3/Vs]
k2 = 3.35; % pump 2 constant [cm^3/Vs]

%% Continuous-Time Linearized Model
% time constants
T1 = (A1/a1) * sqrt(2*h10/g);
T2 = (A2/a2) * sqrt(2*h20/g);
T3 = (A3/a3) * sqrt(2*h30/g);
T4 = (A4/a4) * sqrt(2*h40/g);

Ac = [-1/T1,      0,  A3/(A1*T3),           0;
          0,  -1/T2,           0,  A4/(A2*T4);
          0,      0,       -1/T3,           0;
          0,      0,           0,      -1/T4];

Bc = [    gamma1*k1/A1,                0;
                     0,     gamma2*k2/A2;
                     0, (1-gamma2)*k2/A3;
      (1-gamma1)*k1/A4,               0];

Cc = [1, 0, 0, 0;  % assume the heights in 1 and 2 are directly accessible
      0, 1, 0, 0];

Dc = zeros(2, 2);

%% Discretization
% convert to discrete time using zero-order hold
sys_c = ss(Ac, Bc, Cc, Dc);
sys_d = c2d(sys_c, Ts, 'zoh');
A = sys_d.A;
K = sys_d.B; % K is the input-to-state matrix
C = sys_d.C;

%% Noise Matrices (B and D)
n = 4; % states
p = 2; % outputs

% define arbitrary variances for the physical process and sensors
% joint covariance matrix Sigma

% choice 1: small noise
Sigma = randn(6) * 0.5e-2;
% choice 2: big noise
% Sigma = randn(6) * 1e-1;

% ensure Sigma is symmetric positive definite
Sigma = Sigma * Sigma';

% Cholesky factorization to find B and D
NoiseMat = chol(Sigma, 'lower');
B = NoiseMat(1:n, :);
D = NoiseMat(n+1:end, :);

%% Display and Save
disp('---------------------------------------------------------');
disp(['Matrix [B; D] size: ', num2str(size(NoiseMat, 1)), ' x ', num2str(size(NoiseMat, 2))]);
disp(['Rank of [B; D]: ', num2str(rank(NoiseMat))]);
disp('---------------------------------------------------------');

save('quadruple_tank.mat', 'A', 'K', 'B', 'C', 'D', 'Sigma', 'Ts');
disp('Model saved to quadruple_tank.mat');
end