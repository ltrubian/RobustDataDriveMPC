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
%       mpc:    which strategy to use the MPC
%               RKF-ext: exted the model to the N time horizon and make RKF
%                        to that extended model (just starting point x0 is
%                        given to MPC)
%               RKF:     compute the RKF on the nominal model (just
%                        starting point x0 is given to MPC)
%               LFM:     the time-varying LFM is computed and used for the
%                        prediction x0 (LFM model and x0 are given to MPC)
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
    con_params.beta     (1,1) double {mustBeBetween(con_params.beta,0,1)}
    con_params.mpc      (1,1) string {mustBeMember(con_params.mpc,["RKF", "RKF-ext", "LFM"])} = "RKF"
    con_params.options
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
    err = 0;
    for k=max(t-con_params.L, 1):(t-1) % time-window L
        % forgetting factor beta
        err = err * con_params.beta + norm(simY(k) - model_con.C * filt.x_pred(k),2)^2;
    end
    optimal_values(ff) = optimal_values(ff) + err;
end

for ff = 1:length(filters)
    filt = filters(ff);

    switch con_params.mpc
        case {"RKF", "RKF-ext"}
            % the only difference between the two is how they compute the
            % starting point x0 for the MPC, that is the one-step ahead
            % prediction
            switch con_params.mpc
                case "RKF-ext"
                    % assembly the extended model struct
                    model_ext = struct( ...
                        "A", kron(eye(N), model_con.A), ...
                        "B", kron(eye(N), model_con.B), ...
                        "C", kron([1, zeros(1,N-1)], model_con.C), ...
                        "D", kron([1, zeros(1,N-1)], model_con.D), ...
                        "K", kron(eye(N), model_con.K));
                    % compute the one step ahead prediction for model ext
                    [x0, filters(ff).Vex] = ...
                        RobustKalmanFilter(model_ext, ...
                        filters(ff).Vex, ...
                        paddata(filt.x_pred(:,t),n*N,Side="trailing"), ...
                        simY(:,t), zeros(m*N,1), filt.c);
                case "RKF"
                    % compute the one step ahead prediction for nominal
                    % model
                    x0 = paddata( ...
                        RobustKalmanFilter(model_con, ...
                        filt.V(:,:,t),filt.x_pred(:,t),simY(:,t),zeros(m,1),filt.c), ...
                        n*N,Side="trailing");
            end

            % compute optimal value
            [optimal_u(:,ff), opt_value_tmp] = MPCOptimizer(x0, ...
                model_con.A, model_con.K, model_con.C, model_con.weights, N, reference(:,t:t+N-1), ...
                model_con.x_min, model_con.x_max,model_con.u_min,model_con.u_max, con_params.options);

        case "LFM"
            % compute the time-varying LFM
            [A, B, C, D] = LeastFavorableModel(model_con, filt.V(:,:,t), filt.c, N);
            % assembly the LF model struct
            model_lfm = struct("A", A(:,:,1), "B", B(:,:,1), ...
                "C", C(:,:,1), "D", D(:,:,1), "K", zeros(2*n,m));
            % Kalman gain to compute the prediction error
            G = (filt.V(:,:,t)*model_con.C' + model_con.B*model_con.D') / ...
                (model_con.C*filt.V(:,:,t)*model_con.C' + model_con.D*model_con.D');
            % use the prediction error and the extended model to compute
            % the prediciont fot the MPC
            [x_pred, filters(ff).Vlfm(:,:,t+1)] = RobustKalmanFilter(model_lfm, ...
                filt.Vlfm(:,:,t),[filt.x_pred(:,t); ...
                G * (simY(:,t) - model_con.C*filt.x_pred(:,t))] , ...
                simY(:,t),zeros(m,1),filt.c);

            x0 = sparse(1:2*n, 1, x_pred, 2*n*N, 1);    % paddata
            A = A(:,:,2:end);                           % select the A's
            % compute the optimal input
            [optimal_u(:,ff), opt_value_tmp] = MPCOptimizer(x0, ...
                A, [model_con.K; sparse(n,m)], C, model_con.weights, ...
                N, reference(:,t:t+N-1), ...
                [model_con.x_min; repmat(-Inf,n,1)], [model_con.x_max; repmat(+Inf,n,1)], ...
                model_con.u_min,model_con.u_max, con_params.options);
    end
    % update the optimal values
    optimal_values(ff) = optimal_values(ff) + opt_value_tmp;
end

% select the optimal input
[~, c_index] = min(optimal_values);
simU = optimal_u(:,c_index);
end