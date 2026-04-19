function [simX, simY, simU, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
    steps_sim, init_con, reference, set_c, debug, con_params)
%LOOPSIMULATION Simulate cloosed-loop system
%
%       <usage here>
%
% INPUT:
%   model_sim:  model to simulate
%   model_con:  nominal startgin model used by MPC
%   steps_sim:  number of step to simulate
%   init_con:   initial condition
%   reference:  reference signal
%   set_c:      set of hyperparamter 'c' to choose from
%   debug:      remove all noise leaving deterministic evolution
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
%   simX:       simulated states
%   simY:       simulated output
%   simU:       controlled input
%   cpuT:       cpu time of the controller
%   filters:    struct with the dynamincs of the set of filters
%   c_index:    the sequence of c's selected by the controller

arguments
    model_sim   (1,1) struct
    model_con   (1,1) struct
    steps_sim   (1,1) double {mustBeInteger(steps_sim)}
    init_con    (:,1) double
    reference   (:,:) double
    set_c       (:,1) double
    debug       (1,1) logical
    con_params.N        (1,1) double
    con_params.L        (1,1) double
    con_params.beta     (1,1) double
    con_params.lfm      (1,1) logical
    con_params.steps    (1,1) double
    con_params.mpc      (1,1) string
    con_params.options
end

n = size(model_sim.A,1);
p = size(model_sim.C,1);
m = size(model_sim.K,2);

%% actual state and control trajectories
simX = NaN(n, steps_sim + 1);
simY = NaN(p, steps_sim + 1);
simU = NaN(m, steps_sim);
cpuT = NaN(size(simU,1), 1);

% initialize the simulation initial information
simX(:,1) = init_con;

%% Filters collection of information
filters = repmat( ...
    struct( ...
        "c", 1, ...                              % placeholder for c value
        "P", repmat(eye(n),1,1,steps_sim+1), ... % nominal initial covariance  = I
        "V", repmat(eye(n),1,1,steps_sim+1), ... % perturbed initial covariance= I
        "lambda", zeros(1,steps_sim+1), ...      % lagrange multipliers
        "x_pred", zeros(n,steps_sim+1), ...      % prediction initial estimate
        "Vex", eye(n*con_params.N), ...  % nominal V matrix for big leap estimation
        "xi", zeros(2*n,1), ...
        "Vlfm", eye(2*n) ...
    ), ...
    size(set_c,1),1);
% initialize values of c
for i=1:size(set_c,1)
    filters(i).c = set_c(i);
end
% sequence of the (indexes) tollerances selected by the controller
c_index = ones(steps_sim,1);

% vt:   noise at time t
vt = randn(n+p,1) * (1 - debug);
if debug
    for i=1:size(set_c,1)
        filters(i).x_pred(:,1) = init_con;
    end
end
try
for t = 1:steps_sim
    %% Output of the system
    simY(:,t) = model_sim.C * simX(:,t) ...
        + model_sim.D * vt;
        % to add the input -> output dynamics, make sure matrix and MPC can
        % deal with it. at the moment MPC is not ready
        % ... + model_sim.J * simU(:,t);


    %% Controller and Filter:
    tic;
    % 1) update input of the system
    [simU(:,t), c_index(t), filters] = Controller(model_con, reference, simY, filters, t, ...
        L=con_params.L, N=con_params.N, beta=con_params.beta, ...
        lfm=con_params.lfm, steps=con_params.steps, mpc=con_params.mpc, options=con_params.options);

    % 2) update filters prediction of the next state
    for ff=1:size(set_c,1)
        [filters(ff).x_pred(:,t+1), filters(ff).V(:,:,t+1), ~, filters(ff).P(:,:,t+1), filters(ff).lambda(1,t+1)] = ...
            RobustKalmanFilter(model_con, ...
            filters(c_index(t)).V(:,:,t), ...      % each filter uses the best
            filters(c_index(t)).x_pred(:,t), ...   % prediction till now
            simY(:,t), simU(:,t), filters(ff).c);
    end
    cpuT(t) = toc;

    %% Simulate the system
    % update noise
    vt = randn(n+p,1) * (1 - debug);
    % update state
    simX(:,t+1) = model_sim.A * simX(:,t) ...   % state
        + model_sim.B * vt ...                  % noise
        + model_sim.K * simU(:,t);              % input

    % fprintf('It: %3d/%3d  CPU time: %2.2f TOTAL time: %4.2f  c: %.1e  lambda: %.2e  pred err: %3.3f\n', ...
    %     t, steps_sim, cpuT(t), sum(cpuT), filters(c_index(t)).c, filters(c_index(t)).lambda(t), norm(filters(c_index(t)).x_pred(:,t) - simX(:,t) ,2));
end
catch e
    disp(getReport(e))
    return
end
end