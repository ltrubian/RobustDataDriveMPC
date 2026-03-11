function [x_next, P_next, V_next, lambda] = robust_kalman_step(x_hat, V_t, y_t, u_t, A, B, C, D, c)
    % 1. Calcolo del guadagno (Usa la covarianza distorta V_t precedente)
    % Equazione: G_t = (A*V_t*C' + B*D') / (C*V_t*C' + D*D')
    tmp = (C * V_t * C' + D * D');
    G_t = (A * V_t * C' + B * D') / tmp;
    
    % 2. Aggiornamento dello stato (Struttura nominale, ma basata su V)
    x_next = A * x_hat + B * u_t + G_t * (y_t - C * x_hat);
    
    % 3. Aggiornamento della covarianza nominale (Equazione di Riccati)
    % Prof's formulation   P_next = A * V_t * A' + B * B' - G_t * Ry * G_t';
    P_next = (A-G_t*C) * V_t * (A-G_t*C)' + (B-G_t*D) * (B-G_t*D)';

    % 4. Risoluzione numerica per lambda (Bisezione)
    % Cerchiamo theta tale che gamma(P_next, lambda) = c

    % --- Vincolo di Ammissibilità ---
    % Per garantire che la covarianza robusta V = (P^-1 - theta*I)^-1 sia definita positiva,
    % il parametro di rischio theta deve essere strettamente minore dell'autovalore 
    % minimo di P^-1, ovvero: theta < 1/lambda_max(P). Oltre questo limite, 
    % l'incertezza stimata "esplode", rendendo il gioco minimax non risolvibile 
    % e il filtro numericamente instabile (varianza negativa).
    lambda = LagrangeMultiplier(P_next, c, "fast");
    
    % 5. Calcolo della covarianza distorta V_next
    % V = (P^-1 - theta*I)^-1 
    V_next = inv( inv(P_next)-eye(size(A))/lambda );
end
