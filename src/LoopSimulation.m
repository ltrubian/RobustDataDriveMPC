function [simX, simY, simU, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
    steps_sim, init_con, reference, set_c, con_params)
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
    con_params.N        (1,1) double
    con_params.L        (1,1) double
    con_params.beta     (1,1) double
    con_params.lfm      (1,1) logical
    con_params.steps    (1,1) double
    con_params.mpc      (1,1) string
end

n = size(model_sim.A,1);
p = size(model_sim.C,1);
m = size(model_sim.B,2);

%% actual state and control trajectories
simX = zeros(n, steps_sim + 1);
simY = zeros(p, steps_sim + 1);
simU = zeros(m, steps_sim);
cpuT = zeros(size(simU,1), 1);

% initialize the simulation initial information
simX(:,1) = init_con;

%% Filters collection of information
filters = repmat( ...
    struct( ...
        "c", 1, ...                              % placeholder for c value
        "P", repmat(eye(n),1,1,steps_sim+1), ... % nominal initial covariance  = I
        "V", repmat(eye(n),1,1,steps_sim+1), ... % perturbed initial covariance= I
        "lambda", zeros(1,steps_sim+1), ...      % lagrange multipliers
        "x_pred", zeros(n,steps_sim+1) ...
    ), ...                                       % prediction initial estimate
    size(set_c,1),1);
% initialize values of c
for i=1:size(set_c,1)
    filters(i).c = set_c(i);
end
% sequence of the (indexes) tollerances selected by the controller
c_index = ones(steps_sim,1);

for t = 1:steps_sim
    %% Output of the system
    simY(:,t) = model_sim.C * simX(:,t) + model_sim.D * randn(p,m);

    %% Controller and Filter:
    tic;
    % 1) update input of the system
    [simU(:,t), c_index(t)] = Controller(model_con, reference, simY, filters, t, ...
        L=con_params.L, N=con_params.N, beta=con_params.beta, ...
        lfm=con_params.lfm, steps=con_params.steps, mpc=con_params.mpc);

    % 2) update filters prediction of the next state
    for ff=1:size(set_c,1)
        [filters(ff).x_pred(:,t+1), filters(ff).V(:,:,t+1), ~, filters(ff).P(:,:,t+1), filters(ff).lambda(1,t+1)] = ...
            RobustKalmanFilter(model_con, ...
            filters(ff).V(:,:,t), ...      % each filter uses its own previous
            filters(ff).x_pred(:,t), ...   % prediction: parallel approach
            simY(:,t), simU(:,t), filters(ff).c);
    end
    cpuT(t) = toc;

    %% Simulate the system
    % at the moment the matrix B is the same for the input and the noise.
    % This is NOT in general the case
    simX(:,t+1) = model_sim.A * simX(:,t) + ...
        model_sim.B * simU(:,t) + 0.05 * randn(2,1); % model_sim.B * randn(m,1);

    fprintf('\rIt: %3d/%3d  CPU time: %2.2f TOTAL time: %4.2f  c: %.1e  lambda: %.2e  pred err: %3.3f', ...
        t, steps_sim, cpuT(t), sum(cpuT), filters(c_index(t)).c, filters(c_index(t)).lambda(t), norm(filters(c_index(t)).x_pred(:,t) - simX(:,t) ,2));
end
fprintf("\n");
end