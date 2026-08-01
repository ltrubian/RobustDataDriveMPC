function generate_quadruple_tank_model()
    % GENERATE_QUADRUPLE_TANK_MODEL 
    % Generates the discrete-time state-space matrices A, B, C, D, K
    % for the Quadruple-Tank Process (Johansson, 2000).
    %
    % The resulting model matches the structure:
    % x(t+1) = A x(t) + B v(t) + K u(t)
    % y(t)   = C x(t) + D v(t)
    % where the joint noise covariance [B; D] is square and invertible.
    
    disp('Generating Quadruple-Tank Model...');
    
    %% Configuration
    Ts = 1; % Sampling time in seconds
    
    %% Physical Parameters (Johansson, 2000)
    A1 = 28; A3 = 28; % Area of tanks 1 and 3 [cm^2]
    A2 = 32; A4 = 32; % Area of tanks 2 and 4 [cm^2]
    a1 = 0.071; a3 = 0.071; % Cross-section of outlet hole 1 and 3 [cm^2]
    a2 = 0.057; a4 = 0.057; % Cross-section of outlet hole 2 and 4 [cm^2]
    g = 981; % Gravitational constant [cm/s^2]
    
    %% Operating Point
    gamma1 = 0.70; gamma2 = 0.60;
    h10 = 12.4; h20 = 12.7; h30 = 1.8; h40 = 1.4; % Initial conditions [cm]
    k1 = 3.33; % Pump 1 constant [cm^3/Vs]
    k2 = 3.35; % Pump 2 constant [cm^3/Vs]
    
    %% Continuous-Time Linearized Model
    % Time constants
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
          
    Cc = [1, 0, 0, 0;  % Assume the heights in 1 and 2 are directly accessible
          0, 1, 0, 0];
          
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
    
    % Define arbitrary variances for the physical process and sensors
    % Joint covariance matrix Sigma
    
    % Choice 1: small noise
    Sigma = randn(6) * 0.5e-2;
    % Choice 2: big noise
    % Sigma = randn(6) * 1e-1;
    
    % Ensure Sigma is symmetric positive definite
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