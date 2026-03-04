function [x_next, P_next, V_next, theta] = robust_kalman_step(x_hat, V_t, y_t, u_t, A, B, C, D, c)
    % 1. Calcolo del guadagno (Usa la covarianza distorta V_t precedente)
    % Equazione: G_t = (A*V_t*C' + B*D') / (C*V_t*C' + D*D')
    Ry = (C * V_t * C' + D * D');
    G_t = (A * V_t * C' + B * D') / Ry;
    
    % 2. Aggiornamento dello stato (Struttura nominale, ma basata su V) [cite: 56]
    x_next = A * x_hat + B * u_t + G_t * (y_t - C * x_hat);
    
    % 3. Aggiornamento della covarianza nominale (Equazione di Riccati)
    % P_next = A*V_t*A' + B*B' - G_t * Ry * G_t'
    P_next = A * V_t * A' + B * B' - G_t * Ry * G_t';
    
    % 4. Risoluzione numerica per theta (Bisezione)
    % Cerchiamo theta tale che gamma(P_next, theta) = c
    % Limite superiore: theta < 1/max(eig(P_next)) per garantire definita positività [cite: 59, 62]
    max_eig_P = max(eig(P_next));
    theta_max = 1 / (max_eig_P + 1e-6);
    theta_min = 0;
    
    % Metodo della bisezione
    for i = 1:50
        theta_mid = (theta_min + theta_max) / 2;
        val = calculate_gamma(P_next, theta_mid);
        if val < c
            theta_min = theta_mid;
        else
            theta_max = theta_mid;
        end
    end
    theta = theta_min;
    
    % 5. Calcolo della covarianza distorta V_next
    % V = (P^-1 - theta*I)^-1 
    V_next = inv(inv(P_next) - theta * eye(size(P_next)));
end
