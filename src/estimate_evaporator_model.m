% estimate_evaporator_model.m
% Identifies a state-space model for the evaporator dataset.
% The model structure follows the equation:
% x(t+1) = A x(t) + B v(t) + K u(t)
% y(t)   = C x(t) + D v(t)
% where the joint noise covariance [B; D] is square and invertible.
%
% This uses a deterministic-stochastic Subspace Identification approach
% to ensure the resulting process/measurement noise covariance is full rank.

clear; close all; clc;

%% 1. Load Data
% Check if the uncompressed file exists, otherwise uncompress it
if ~isfile('evaporator.dat')
    if isfile('evaporator.dat.gz')
        gunzip('evaporator.dat.gz');
    else
        error('Data file evaporator.dat or evaporator.dat.gz not found.');
    end
end

disp('Loading dataset...');
data = load('evaporator.dat');

% Center the data around the operating point (remove means)
u = data(:, 1:3);
y = data(:, 4:6);
u = bsxfun(@minus, u, mean(u));
y = bsxfun(@minus, y, mean(y));

N = size(u, 1);
m = size(u, 2); % Inputs (3)
p = size(y, 2); % Outputs (3)

%% 2. User-Defined Subspace Parameters
n = 4;   % State dimension
f = 10;  % Past and future horizon lengths (tuning parameter)

disp(['Identifying model with n=', num2str(n), ', m=', num2str(m), ', p=', num2str(p)]);

%% 3. Form Block Hankel Matrices
disp('Forming block Hankel matrices...');
U = zeros(2*f*m, N - 2*f + 1);
Y = zeros(2*f*p, N - 2*f + 1);

for i = 1:2*f
    U((i-1)*m+1:i*m, :) = u(i:N-2*f+i, :)';
    Y((i-1)*p+1:i*p, :) = y(i:N-2*f+i, :)';
end

% Partition into past and future
Up = U(1:f*m, :);
Uf = U(f*m+1:end, :);
Yp = Y(1:f*p, :);
Yf = Y(f*p+1:end, :);

% Instrument variable (past inputs and outputs)
Wp = [Up; Yp];

%% 4. Oblique Projection and SVD (Estimate State Sequence)
disp('Projecting and estimating state sequence (N4SID)...');

% We want to project Yf along Uf onto Wp.
% Mathematically: Yf /_{Uf} Wp. We use least squares.
% Regress Yf onto [Wp; Uf]
Regressors = [Wp; Uf];
Theta = Yf * pinv(Regressors);

% Extract the part corresponding to Wp
L_W = Theta(:, 1:size(Wp,1));

% Compute the projection
Proj = L_W * Wp;

% Perform Singular Value Decomposition on the projection
[U_svd, S_svd, V_svd] = svd(Proj, 'econ');

% Truncate to state dimension 'n'
S_n = S_svd(1:n, 1:n);
V_n = V_svd(:, 1:n);

% Compute the state sequence for the future horizon
% Size: n x (N - 2*f + 1)
X = sqrt(S_n) * V_n'; 

%% 5. Estimate System Matrices via Least Squares
disp('Estimating system matrices A, K, C...');

% States at time t and t+1
X_t   = X(:, 1:end-1);
X_tp1 = X(:, 2:end);

% Inputs and outputs at time t (which corresponds to the first block of Uf and Yf)
u_t = Uf(1:m, 1:end-1);
y_t = Yf(1:p, 1:end-1);

% Solve for A and K: X_{t+1} = A*X_t + K*u_t + w_t
AK = X_tp1 * pinv([X_t; u_t]);
A = AK(:, 1:n);
K = AK(:, n+1:end);

% Solve for C: y_t = C*X_t + e_t
C = y_t * pinv(X_t);

%% 6. Compute Residuals and Noise Matrices
disp('Computing full-rank joint noise covariance and matrices B, D...');

% Compute residuals
w_t = X_tp1 - (A * X_t + K * u_t);
e_t = y_t - C * X_t;

% Joint sample covariance matrix
% Note: Since we derived states from noisy data, w_t and e_t are not perfectly 
% collinear. Thus, Sigma will generically have full rank (n+p).
Sigma = ([w_t; e_t] * [w_t; e_t]') / size(w_t, 2);

% Check condition number to ensure it's numerically invertible
if cond(Sigma) > 1e12
    warning('The sample covariance matrix Sigma is ill-conditioned. Adding a tiny regularization term.');
    Sigma = Sigma + eye(size(Sigma)) * 1e-8;
end

% Perform Cholesky factorization: Sigma = L * L'
% L will be a square lower triangular matrix of size (n+p) x (n+p)
NoiseMat = chol(Sigma, 'lower');

% Extract B and D matrices
B = NoiseMat(1:n, :);
D = NoiseMat(n+1:end, :);

disp('---------------------------------------------------------');
disp(['Matrix [B; D] size: ', num2str(size(NoiseMat, 1)), ' x ', num2str(size(NoiseMat, 2))]);
disp(['Rank of [B; D]: ', num2str(rank(NoiseMat))]);
disp('---------------------------------------------------------');

%% 7. Save the Resulting Matrices
save('evaporator_model.mat', 'A', 'K', 'B', 'C', 'D', 'Sigma');
disp('Model successfully saved to evaporator_model.mat.');
