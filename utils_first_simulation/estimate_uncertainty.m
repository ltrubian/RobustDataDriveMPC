function c_best = estimate_uncertainty(y_hist, u_hist, A, B, C, D, c_candidates, decay)
    num_c = length(c_candidates);
    errors = zeros(1, num_c);
    window = min(length(y_hist)-1, 10); % Finestra temporale
    
    for i = 1:num_c
        c_test = c_candidates(i);
        % Reset filtri locali per il test
        x_tmp = [0;0]; V_tmp = eye(2)*0.1;
        total_err = 0;
        
        for k = (length(y_hist)-window):(length(y_hist)-1)
            % Predizione a un passo
            y_pred = C * x_tmp;
            err = (y_hist(k) - y_pred)^2;
            total_err = total_err * decay + err;
            
            % Update filtro (RKF semplificato per test rapido)
            [x_tmp, ~, V_tmp, ~] = robust_kalman_step(x_tmp, V_tmp, y_hist(k), u_hist(k), A, B, C, D, c_test);
        end
        errors(i) = total_err;
    end
    [~, idx] = min(errors);
    c_best = c_candidates(idx);
end