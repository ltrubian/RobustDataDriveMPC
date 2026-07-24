function [simX, simY, trueY, simU, cpuT, RKFs, c_best] = LoopSimulation(model_sim, model_nom, ...
    steps_sim, init_con, reference, set_c, debug, verbose, con_params)
%LOOPSIMULATION Simulate cloosed-loop system
%
%       <usage here>
%
% INPUT:
%   model_sim:  model to simulate
%   model_nom:  nominal starting model used by MPC
%   steps_sim:  number of step to simulate
%   init_con:   initial condition
%   reference:  reference signal
%   set_c:      set of hyperparamter 'c' to choose from
%   debug:      remove all noise leaving deterministic evolution
%   con_params: see NAMED-VALUE INPUTS of Controller function
%
% OUTPUT:
%   simX:       simulated states
%   simY:       simulated output
%   simU:       controlled input
%   cpuT:       cpu time of the controller
%   filters:    struct with the dynamincs of the set of filters
%   c_index:    the sequence of c's selected by the controller

arguments
    model_sim   (1,1) struct
    model_nom   (1,1) struct
    steps_sim   (1,1) double {mustBeInteger(steps_sim)}
    init_con    (:,1) double
    reference   (:,:) double
    set_c       (:,1) double
    debug       (1,1) logical
    verbose     (1,1) logical
    con_params  (1,1) struct
end

n = size(model_sim.A,1);        % state real world
p = size(model_sim.C,1);        % input real world
m = size(model_sim.K,2);        % output real world
r = size(model_nom.A,1) - n;    % fictitious disturbances (if introduced)
N = con_params.N;

%% actual state and control trajectories
simX = NaN(n, steps_sim + 1);
simY = NaN(p, steps_sim + 1);
simU = NaN(m, steps_sim);
trueY = NaN(p, steps_sim + 1);
cpuT = NaN(size(simU,1), 1);

% initialize the simulation initial information
simX(:,1) = init_con;

%% Filters collection of information
% each type of controller (since it is based to different filter) works 
% with a specific state, state dynamics, and initial covariance
switch con_params.mpc
    case "RKF"
        nc = n+r;       % controll state: state + disturbances (if present)
        nf = n+r;       % estimated state: RKF is the controll sttate
        V_0 = eye(nf);
        model_fil = model_nom;
        model_con = model_nom;
    case "RKF-ext"
        nf = n*N;
        V_0 = kron(ones(N)+eye(N)/10,eye(n));
        model_fil = struct( ...
            "A", kron(eye(N), model_nom.A), ...
            "B", kron(eye(N), model_nom.B), ...
            "C", kron([1, zeros(1,N-1)], model_nom.C), ...
            "D", kron([1, zeros(1,N-1)], model_nom.D), ...
            "K", cell2mat(arrayfun(@(k) (model_nom.A)^k, 0:N-1, 'UniformOutput', false)')* model_nom.K);
        model_con = model_nom;
    case "LFM"
        nc = n+r;       % controll state: state + disturbances (if present)
        nf = 2*(n+r);   % estimated state: LFM is double of controll sttate
        V_0 = eye(nf);
        % in this case the models used by controller and the filter are 
        % updated at each iteration and it is unnecessary to initialize
        % them
        model_fil = NaN;
        model_con = NaN;
end
RKFs = repmat( ...
    struct( ...
    "c", 1, ...                           % placeholder for c value
    "V", repmat(V_0,1,1,steps_sim+1), ... % perturbed initial covariance
    "x_pred", zeros(nf,steps_sim+1)), ... % prediction initial estimate
    size(set_c,1),1);
% initialize values of c
for i=1:size(set_c,1)
    RKFs(i).c = set_c(i);
end
% sequence of the (indexes) tollerances selected by the controller
c_best = ones(steps_sim+1,1);
x_best = RKFs(1).x_pred(:,1);
V_best = RKFs(1).V(:,:,1);

% vt:   noise at time t
vt = randn(n+p,1) * (1 - debug);
if debug
    for i=1:size(set_c,1)
        RKFs(i).x_pred(:,1) = paddata(init_con, nf, Side="trailing");
    end
end
for t = 1:steps_sim
    %% Output of the system
    trueY(:,t) = model_sim.C * simX(:,t);
    simY(:,t) = model_sim.C * simX(:,t) ...
        + model_sim.D * vt;
    % to add the input -> output dynamics, make sure matrix and MPC can
    % deal with it. at the moment MPC is not ready
    % ... + model_sim.J * simU(:,t);


    %% Controller and Filter:
    tic;
    optimal_values = zeros(size(set_c,1),1);
    optimal_u = zeros(N*m, size(set_c,1));

    for cj=1:length(RKFs)
        % the proposed controller need to update both the model used by the
        % filter for the one-step ahead prediction of the free evolution,
        % and the model used by the MPC. This does NOT happen in the other
        % cases since they keep the nominal model for the MPC (and a static
        % model for the filter as RKF-ext)
        if strcmp("LFM", con_params.mpc)
            [A, B, C, D] = LeastFavorableModel(model_nom, RKFs(c_best(t)).V(1:nc,1:nc,t), RKFs(cj).c, N);
            model_fil = struct( ...
                "A", A(:,:,1), "B", B(:,:,1), ...
                "C", C(:,:,1), "D", D(:,:,1), "K", [model_nom.K; zeros(nc,m)]);
            model_con = struct( ...
                "A", A(:,:,2:end), "C", C, "K", [model_nom.K; sparse(nc,m)], ...
                "x_min", [model_nom.x_min; -Inf(nc,1)], ...
                "x_max", [model_nom.x_max; +Inf(nc,1)], ...
                "u_min", model_nom.u_min, "u_max", model_nom.u_max);
        end
        % 1) Prediction step
        [RKFs(cj).x_pred(:,t+1), RKFs(cj).V(:,:,t+1)] = ...
            RobustKalmanFilter(model_fil, ...
            RKFs(c_best(t)).V(:,:,t), ...   % each filter uses the best
            RKFs(c_best(t)).x_pred(:,t), ...   % prediction till now
            simY(:,t), RKFs(cj).c);

        % 2) Controller step
        x0 = paddata(RKFs(cj).x_pred(:,t+1), nf*N, Side="trailing");
        [optimal_u(:,cj), optimal_values(cj)] = MPCOptimizer(x0, ...
            model_con.A, model_con.K, model_con.C, model_nom.weights, ...
            N, reshape(reference(:,t:t+N-1),[],1), ...
            model_con.x_min, model_con.x_max, ...
            model_con.u_min, model_con.u_max, con_params.options);
    end
    optimal_values = optimal_values ./ (eps+max(optimal_values));
    optimal_values_ff = zeros(size(set_c,1),1);
    for cj=1:length(RKFs)
        % 3) Add past prediction error
        err = 0;
        for k=max(t-con_params.L, 1):(t-1) % time-window L
            % forgetting factor beta
            err = err * con_params.beta ...
                + norm(simY(:,k) - model_nom.C * RKFs(cj).x_pred(1:nc,k),2)^2;
        end
        optimal_values_ff(cj) = err;
    end
    optimal_values = optimal_values + (optimal_values_ff ./(eps+max(optimal_values_ff)) );
    % 4) Optimization step
    [~, c_best(t+1)] = min(optimal_values);
    % Update prediction of each filter based on the selected output
    for cj=1:length(RKFs)
        RKFs(cj).x_pred(:,t+1) = RKFs(cj).x_pred(:,t+1) ...
            + model_fil.K * optimal_u(1:m,c_best(t+1));
    end
    % return input
    simU(:,t) = optimal_u(1:m, c_best(t+1));
    cpuT(t) = toc;

    %% Simulate the system
    % update noise
    vt = randn(n+p,1) * (1 - debug);
    % update state
    simX(:,t+1) = model_sim.A * simX(:,t) ...   % state
        + model_sim.B * vt ...                  % noise
        + model_sim.K * simU(:,t);              % input

    if verbose
        fprintf(['It: %3d/%3d  CPU: %2.2f TOTAL: %4.2f  c: %.1e ' ...
            'pred err: %3.3f min_opt: %3.3f max_opt: %3.3f\n'], ...
            t, steps_sim, cpuT(t), sum(cpuT), RKFs(c_best(t)).c, ...
            norm(RKFs(c_best(t)).x_pred(1:n,t) - simX(:,t) ,2), ...
            min(optimal_values),max(optimal_values));
    end
end
end