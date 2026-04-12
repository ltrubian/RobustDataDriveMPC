function [simU, c_index, filters] = Controller(model_con, reference, simY, filters, t, con_params)
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
%       lfm:    apply Least-Favorable Model (true/false)
%       steps:  1 to combine estimation and controller;
%               2 to have estimation minimizing before and then controller
%       mpc:    which matlab function to use for the MPC controller
%               fmincon: more readble but slow (not suitable for big N)
%               quadprog: fast quadratic solver for sparse mpc
%                         implemntation
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
    con_params.N        (1,1) double {mustBeInteger(con_params.N)}
    con_params.L        (1,1) double {mustBeInteger(con_params.L)}
    con_params.beta     (1,1) double %{mustBeBetween(con_params.beta,0,1)}
    con_params.lfm      (1,1) logical = true
    con_params.steps    (1,1) double {mustBeMember(con_params.steps,[1,2])} = 1
    con_params.mpc      (1,1) string {mustBeMember(con_params.mpc,["fmincon","quadprog", "quad-pro"])} = "quadprog"
    con_params.options  (1,1) optim.options.Quadprog
end

n_filts = length(filters);
n = size(model_con.A,1);
m = size(model_con.K,2);
N = con_params.N;
assert(m == 1, "the call fmincon for vectorial input is NOT yet ready")

% store optimal results
optimal_values = zeros(n_filts,1);
optimal_u = zeros(m, n_filts);

% compute uncertainty for each filter
for ff=1:length(filters)
    filt = filters(ff);
    optimal_values(ff) = optimal_values(ff) + past_prediction_error( ...
        simY, filt, model_con.C, t, con_params.L, con_params.beta);
end

% select if use one estimation or leace the choice to the controller cost
% function: steps
%   1 -> one minimization is done summing cost values of estimators and
%        controllers
%   2 -> THE best estimator is choosen to pass to the controller (the
%        following "for-loop" collapses to one iteration)
switch con_params.steps
    case 1
        chosen_filters = 1:length(filters);
    case 2
        [~, c_index] = min(optimal_values);
        chosen_filters = [c_index];
end

options = optimoptions('fmincon', 'Display', 'off');


for ff = chosen_filters
    filt = filters(ff);
    A = model_con.A; B = model_con.K; C = model_con.C;

    switch con_params.mpc
        case "fmincon"
            u_init = zeros(con_params.N, 1);
            cost_func = @(u_seq) mpc_cost(u_seq, filt.x_pred(:,t), ...
                A, B, C, con_params.N, ...
                reference(t:min(t+con_params.N-1, end)));

            % WARNING: non ho la minima idea di come rendere i vincoli di questa
            % funzione validi per input u che siano vettori. OPS
            [u_tmp, opt_value_tmp ] = fmincon(cost_func, u_init, [], [], [], [], ...
                model_con.u_min*ones(con_params.N,1), model_con.u_max*ones(con_params.N,1), [], options);
            % e proprio per questo non so come modificare la seguente riga per
            % input vettoriali. HELP
            optimal_u(:,ff) = u_tmp(1);

        case "quadprog"
            if con_params.lfm
                model_extended.A = kron(eye(N), model_con.A);
                model_extended.B = kron(eye(N), model_con.B);
                model_extended.C = kron([1, zeros(1,N-1)], model_con.C);
                model_extended.D = kron([1, zeros(1,N-1)], model_con.D);
                model_extended.K = kron(eye(N), model_con.K);
                % model_extended.J =
                [x0, filters(ff).Vex] = ...
                    RobustKalmanFilter(model_extended, ...
                    filters(ff).Vex, ...
                    paddata(filt.x_pred(:,t),n*N,Side="trailing"), ...
                    simY(:,t), zeros(m*N,1), filt.c);
            else
                % x0 = paddata(A*filt.x_pred(:,t) ,n*N,Side="trailing");
                x0 = paddata( ...
                    RobustKalmanFilter(model_con, ...
                    filt.V(:,:,t),filt.x_pred(:,t),simY(:,t),zeros(m,1),filt.c), ...
                    n*N,Side="trailing");
            end

            [optimal_u(:,ff), opt_value_tmp] = MPCOptimizer(x0, ...
                A, B, C, model_con.weights, N, reference(:,t:t+N-1), ...
                model_con.x_min, model_con.x_max,model_con.u_min,model_con.u_max, con_params.options);
        case "quad-pro"
            [A, ~, C] = LeastFavorableModel(model_con, filt.V(:,:,t), filt.c, N);
            x0 = paddata( ...
                    RobustKalmanFilter(model_con, ...
                    filt.V(:,:,t),filt.x_pred(:,t),simY(:,t),zeros(m,1),filt.c), ...
                    2*n*N,Side="trailing");
            A = A(:,:,1:end-1);
            [optimal_u(:,ff), opt_value_tmp] = MPCTV(x0, ...
                A, [B; sparse(n,m)], C, model_con.weights, N, reference(:,t:t+N-1), ...
                ... [model_con.x_min; repmat(-inf,n,1)], [model_con.x_max; repmat(+inf,n,1)], ...
                model_con.x_min, model_con.x_max, ...
                model_con.u_min,model_con.u_max, con_params.options);
    end
    optimal_values(ff) = optimal_values(ff) + opt_value_tmp;
end

switch con_params.steps
    case 1
        [~, c_index] = min(optimal_values);
    case 2
        % notice that the c_index was already selected before
end

simU = optimal_u(:,c_index);

end


function J = mpc_cost(u_seq, x0, A, B, C, N, y_ref)
J = 0;
x_k = x0;

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