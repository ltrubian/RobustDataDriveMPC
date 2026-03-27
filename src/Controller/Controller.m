function [simU, c_index] = Controller(model_con, reference, simY, filters, t, con_params)
%CONTROLLER compute the inputs
%
%       <usage here>
%
% INPUT:
%   model_con:  nominal model used by MPC
%   reference:  reference signal
%   simY:       vector of observations of y
%   filters:    vector of available filters
%   t:          current time instant (index of different vectors)
%
% NAMED-VALUE INPUTS:
%   con_params:
%       N:      prediction horizon of MPC
%       L:      time windows relevant for estimation
%       beta:   forgetting factor
%
% OUTPUT:
%   simU:       input of the system
%   c_index:    index of the best value of c according the whole controller

arguments
    model_con   (1,1) struct
    reference   (:,:) double
    simY        (:,:) double
    filters     (:,:) struct
    t           (1,1) {mustBeInteger(t)}
    con_params.N    (1,1) {mustBeInteger(con_params.N)}
    con_params.L    (1,1) {mustBeInteger(con_params.L)}
    con_params.beta (1,1) double
end

n_filts = length(filters);
m = size(model_con.B,2);
assert(m == 1, "the call fmincon for vectorial input is NOT yet ready")

options = optimoptions('fmincon', 'Display', 'off');

% store optimal results
optimal_values = zeros(n_filts,1);
optimal_u = zeros(m, n_filts);

for ff=1:n_filts
    filt = filters(ff);

    %% MPC controller step fmincon implementation
    % %%%%%%%%%%% start %%%%%%%%%%%
    % u_init = zeros(con_params.N, 1);
    % cost_func = @(u_seq) mpc_cost(u_seq, filt.x_pred(:,t), filt.P(:,:,t), ...
    %     model_con.A, model_con.B, model_con.C, con_params.N, ...
    %     reference(t:min(t+con_params.N-1, end)), filt.lambda(t));
    %
    % % WARNING: non ho la minima idea di come rendere i vincoli di questa
    % % funzione validi per input u che siano vettori. OPS
    % [u_tmp, optimal_values(ff) ] = fmincon(cost_func, u_init, [], [], [], [], ...
    %     model_con.u_min*ones(con_params.N,1), model_con.u_max*ones(con_params.N,1), [], options);
    % % e proprio per questo non so come modificare la seguente riga per
    % % input vettoriali. HELP
    % optimal_u(:,ff) = u_tmp(1);
    % %%%%%%%%%%% end %%%%%%%%%%%
    %% MPC quadprog implementation
    % %%%%%%%%%%% start %%%%%%%%%%%
    A = model_con.A; C = model_con.C;
    if not(isapprox(filt.lambda(t),0))
        distortion = (eye(size(model_con.A)) - filt.P(:,:,t)/filt.lambda(t));
        A = distortion \ model_con.A;
        C = model_con.C / distortion;
    end
    [optimal_u(:,ff), optimal_values(ff)] = MPCOptimizer(filt.x_pred(:,t), ...
        A,model_con.B,C, model_con.weights, con_params.N, reference(:,t:t+con_params.N-1), ...
        model_con.x_min, model_con.x_max,model_con.u_min,model_con.u_max);
    % %%%%%%%%%%% end %%%%%%%%%%%
    %% uncertainty evaluation for each filter
    % computation of the following kind at the end
    optimal_values(ff) = optimal_values(ff) + past_prediction_error( ...
        simY, filt, model_con.C, t, con_params.L, con_params.beta);
end

% minimization of the combined measures
[~, c_index] = min(optimal_values);
simU = optimal_u(:,c_index);

end


function J = mpc_cost(u_seq, x0, P, A, B, C, N, y_ref, lambda)
J = 0;
x_k = x0;

% Least Favorable Model construction
if not(isapprox(lambda,0))
    distortion = (eye(size(A)) - P/lambda);
    A = distortion \ A;
    C = C / distortion;
end

for k = 1:N
    % Prediction
    y_k = C * x_k;
    ref = y_ref(min(k, end));

    % Costo: Tracking + Energia ingresso
    J = J + norm(y_k - ref)^2 + 0.1 * u_seq(k)^2;

    % Evoluzione (Approssimazione Least Favorable: si potrebbe aggiungere bias theta*P*x)
    % Per semplicità, usiamo la dinamica nominale partendo dallo stato "robustificato"
    x_k = A * x_k + B * u_seq(k);
end
end

function err = past_prediction_error(y, filter, C, t, L, beta)
% with this function the filter uses the prediction it has done at that
% time. that prediction was already computed starting from its own prediction
% parallel approach
err = 0;
for k=max(t-L, 1):(t-1)
    err = err * beta + norm(y(k) - C * filter.x_pred(k),2)^2;
end
end