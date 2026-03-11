%% Framework Robust MPC basato su Robust Kalman Filter (RKF)
clear; clc; close all;
rng(0);
addpath("Controller/")
addpath("RobustKalmanFilter/")
% --- 1. Inizializzazione del Modello ---
% Modello Nominale (es. sistema 2D semplificato)
A = [1.1 1; 0 1]; 
B = [0.5; 1]; 
C = [1 0]; 
D_noise = [0.1];% 0.05]; % Matrice D del paper: lega il rumore v_t a x_t+1 e y_t [cite: 63]

% Modello Reale (con incertezza parametrica)
A_real = [1.15 1; 0.01 1]; 
B_real = [0.48; 1.05];
C_real = C;

% Parametri di Simulazione
T = 300;           % Tempo totale
N = 5;            % Orizzonte MPC
c_candidates = logspace(-6, -3, 10); % Set di possibili valori per c
decay_rate = 0.9; % Fattore di decadimento per stima c

% Stato iniziale
x_real = [1; 0];
x_hat = [0; 0];
P = eye(2) * 0.1; % Covarianza iniziale nominale
V = P;            % Covarianza iniziale distorta (V_0 = P_0) 

% Vincoli MPC
u_min = -2; u_max = 2;
y_ref = ones(1, T) * 5;

% Storico per stima c
y_history = [];
u_history = [];
c_hat = c_candidates(1); % Valore iniziale

% --- 5. Loop Principale di Simulazione ---
history_x = zeros(2, T);
history_u = zeros(1, T);
history_x_hat = zeros(2, T);
history_y_hat = zeros(1, T);

for t = 1:T
    tic;
    % A. Ricezione uscita dal sistema reale
    noise = randn(1, 1);
    y_t = C_real * x_real + 0.1 * noise;
    y_history = [y_history, y_t];
    
    % B. Stima dell'incertezza c (Punto 3 dei requisiti)
    if t > 5
        c_hat = estimate_uncertainty(y_history, u_history, A, B, C, D_noise, c_candidates, decay_rate);
    end
    
    % C. Ottimizzazione MPC (Punto 4 dei requisiti)
    % Nota: Usiamo il "least favorable model" per le predizioni future
    u_opt = mpc_optimizer(x_hat, V, A, B, C, N, y_ref(t:min(t+N-1, T)), u_min, u_max, c_hat);
    u_t = u_opt(1);
    u_history = [u_history, u_t];
    
    % D. Update dello stato reale (Simulazione)
    x_real = A_real * x_real + B_real * u_t + 0.05 * randn(2,1);
    
    % E. Esecuzione Robust Kalman Filter (Punto 2 dei requisiti)
    [x_hat, P, V, theta_t] = robust_kalman_step(x_hat, V, y_t, u_t, A, B, C, D_noise, c_hat);
    
    % Salvataggio della stima dello stato
    history_x_hat(:, t) = x_hat;

    % Salvataggio stima uscita
    history_y_hat(t) = C * x_hat; 

    % Salvataggio dati
    history_x(:, t) = x_real;
    history_u(t) = u_t;
    fprintf('Step %d: c_stimatol = %.3f, theta = %.3f\n', t, c_hat, theta_t);
end

% Plot Risultati
figure;

subplot(3,1,1); 
plot(history_x(1,:)); 
hold on; 
plot(history_x_hat(1,:)); 
legend('True State','RKF Estimate');
title('True VS Estimated States');

subplot(3,1,2); 
plot(y_history(1,:));
hold on;
plot(history_y_hat(1,:));
legend('True Output','RKF Estimate');
title('True VS Estimated Output');

subplot(3,1,3); 
stairs(history_u); 
title('Ingresso di Controllo u');












