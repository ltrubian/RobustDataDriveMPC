function J = mpc_cost(u_seq, x0, V0, A, B, C, N, y_ref, c)
    J = 0;
    x_k = x0;
    V_k = V0;
    % Placeholder per matrici LFM (Least Favorable Model) [cite: 5, 29]
    % In un gioco minimax, il modello peggiore sposta la media o aumenta il rumore.
    % Qui usiamo la stima robusta aggiornata per proiettare il peggior caso.
    
    for k = 1:N
        % Predizione nominale (la robustezza è nel valore di x_k stimato dal RKF)
        y_k = C * x_k;
        ref = y_ref(min(k, length(y_ref)));
        
        % Costo: Tracking + Energia ingresso
        J = J + norm(y_k - ref)^2 + 0.1 * u_seq(k)^2;
        
        % Evoluzione (Approssimazione Least Favorable: si potrebbe aggiungere bias theta*P*x)
        % Per semplicità, usiamo la dinamica nominale partendo dallo stato "robustificato"
        x_k = A * x_k + B * u_seq(k);
    end
end