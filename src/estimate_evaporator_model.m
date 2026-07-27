% estimate_evaporator_model.m
% Identifies a state-space model for the evaporator dataset.
% The model structure follows the equation:
% x(t+1) = A x(t) + B v(t) + K u(t)
% y(t)   = C x(t) + D v(t)
% where the joint noise covariance [B; D] is square and invertible.
%
% This script performs the following steps:
% 1. Divides the dataset (50/50).
% 2. Estimates the model using the first half.
% 3. Validates the model on the second half.
% 4. Re-estimates the model using the full dataset.
% 5. Saves the final model.

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

N_total = size(u, 1);
m = size(u, 2); % Inputs (3)
p = size(y, 2); % Outputs (3)

%% 2. Split Dataset
split_idx = floor(N_total / 2);
u_est = u(1:split_idx, :);
y_est = y(1:split_idx, :);
u_val = u(split_idx+1:end, :);
y_val = y(split_idx+1:end, :);

%% 3. User-Defined Subspace Parameters
n = 4;   % State dimension
f = 10;  % Past and future horizon lengths (tuning parameter)

disp(['Identifying model with n=', num2str(n), ', m=', num2str(m), ', p=', num2str(p)]);

%% 4. Estimate on First Half
disp('---------------------------------------------------------');
disp('Estimating model on the FIRST HALF of the dataset...');
[A_est, K_est, B_est, C_est, D_est, Sigma_est] = identify_model(u_est, y_est, n, f);

%% 5. Validate on Second Half
disp('---------------------------------------------------------');
disp('Validating the estimated model on the SECOND HALF of the dataset...');

N_val = size(u_val, 1);
x_sim = zeros(n, N_val);
y_sim = zeros(p, N_val);

% Initialize state (assume zero since data is mean-centered)
x_sim(:, 1) = zeros(n, 1);

for t = 1:N_val-1
    y_sim(:, t) = C_est * x_sim(:, t);
    x_sim(:, t+1) = A_est * x_sim(:, t) + K_est * u_val(t, :)';
end
y_sim(:, N_val) = C_est * x_sim(:, N_val);
y_sim = y_sim'; % Convert to N_val x p

% Calculate fit percentages
fits = zeros(1, p);
for i = 1:p
    y_true = y_val(:, i);
    y_hat = y_sim(:, i);
    fits(i) = 100 * (1 - norm(y_true - y_hat) / norm(y_true - mean(y_true)));
    disp(['Output y', num2str(i), ' Fit: ', num2str(fits(i), '%.2f'), '%']);
end

% Plotting
figure('Name', 'Validation Results', 'NumberTitle', 'off');
for i = 1:p
    subplot(p, 1, i);
    plot(1:N_val, y_val(:, i), 'k', 1:N_val, y_sim(:, i), 'r--');
    legend('True', 'Simulated');
    title(['Output y', num2str(i), ' (Fit: ', num2str(fits(i), '%.2f'), '%)']);
    xlabel('Time step');
    ylabel(['y_', num2str(i)]);
end

%% 6. Re-estimate on Whole Dataset
disp('---------------------------------------------------------');
disp('Re-estimating model on the WHOLE dataset...');
[A, K, B, C, D, Sigma] = identify_model(u, y, n, f);

%% 7. Save the Resulting Matrices
save('evaporator_model.mat', 'A', 'K', 'B', 'C', 'D', 'Sigma');
disp('---------------------------------------------------------');
disp('Final model successfully saved to evaporator_model.mat.');


%% ========================================================================
%  Helper Function for Subspace Identification (N4SID)
%  ========================================================================
function [A, K, B, C, D, Sigma] = identify_model(u_data, y_data, n, f)
N = size(u_data, 1);
m = size(u_data, 2);
p = size(y_data, 2);

% Form block Hankel matrices
U = zeros(2*f*m, N - 2*f + 1);
Y = zeros(2*f*p, N - 2*f + 1);
for i = 1:2*f
    U((i-1)*m+1:i*m, :) = u_data(i:N-2*f+i, :)';
    Y((i-1)*p+1:i*p, :) = y_data(i:N-2*f+i, :)';
end

Up = U(1:f*m, :);
Uf = U(f*m+1:end, :);
Yp = Y(1:f*p, :);
Yf = Y(f*p+1:end, :);
Wp = [Up; Yp];

% Oblique Projection
Regressors = [Wp; Uf];
Theta = Yf * pinv(Regressors);
L_W = Theta(:, 1:size(Wp,1));
Proj = L_W * Wp;

% SVD to estimate state sequence
[U_svd, S_svd, V_svd] = svd(Proj, 'econ');
S_n = S_svd(1:n, 1:n);
V_n = V_svd(:, 1:n);
X = sqrt(S_n) * V_n';

% Least Squares for A, K, C
X_t   = X(:, 1:end-1);
X_tp1 = X(:, 2:end);
u_t = Uf(1:m, 1:end-1);
y_t = Yf(1:p, 1:end-1);

AK = X_tp1 * pinv([X_t; u_t]);
A = AK(:, 1:n);
K = AK(:, n+1:end);

C = y_t * pinv(X_t);

% Compute Residuals and Noise Matrices
w_t = X_tp1 - (A * X_t + K * u_t);
e_t = y_t - C * X_t;

% Joint sample covariance
Sigma = ([w_t; e_t] * [w_t; e_t]') / size(w_t, 2);

% Check condition number and regularize if necessary
if cond(Sigma) > 1e12
    Sigma = Sigma + eye(size(Sigma)) * 1e-8;
end

% Extract B and D matrices
NoiseMat = chol(Sigma, 'lower');
B = NoiseMat(1:n, :);
D = NoiseMat(n+1:end, :);

disp(['  -> Matrix [B; D] size: ', num2str(size(NoiseMat, 1)), ' x ', num2str(size(NoiseMat, 2))]);
disp(['  -> Rank of [B; D]: ', num2str(rank(NoiseMat))]);
end
