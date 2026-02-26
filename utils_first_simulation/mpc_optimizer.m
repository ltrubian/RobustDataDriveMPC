function u_opt = mpc_optimizer(x0, V0, A, B, C, N, y_ref, u_min, u_max, c)
    u_init = zeros(N, 1);
    options = optimoptions('fmincon', 'Display', 'off');
    
    % Definizione funzione di costo
    cost_func = @(u_seq) mpc_cost(u_seq, x0, V0, A, B, C, N, y_ref, c);
    
    % Ottimizzazione con vincoli
    u_opt = fmincon(cost_func, u_init, [], [], [], [], u_min*ones(N,1), u_max*ones(N,1), [], options);
end