function [simX, simU, cpuT, filters] = LoopSimulation(model_sim, model_con, ...
    steps_sim, ini_con, reference, set_c, con_params)
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
%
% OUTPUT:
%   simX:       simulated states
%   simU:       controlled input
%   cpuT:       cpu time of the controller
%   filters:    struct with the dynamincs of the set of filters
%
% DEV-STATUS:
%   NEVER RUN

arguments
    model_sim   (1,1) struct
    model_con   (1,1) struct
    steps_sim   (1,1) double {mustBeInteger(steps_sim)}
    ini_con     (:,1) double
    reference   (:,:) double
    set_c       (1,:) double
    con_params.N    (1,1) double = 20
    con_params.L    (1,1) double = 10
    con_params.beta (1,1) double = 0.999
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
simX(:,1) = ini_con;

%% Filters collection of information
filters = repmat(struct("c", 1, ...
    "P", repmat(eye(n),1,1,steps_sim+1), ... % nominal initial covariance  = I
    "V", repmat(eye(n),1,1,steps_sim+1), ... % perturbed initial covariance= I
    "lambda", ones(1,steps_sim+1), ...       % lagrange multipliers
    "x_pred", zeros(n,steps_sim+1)), ...     % prediction initial estimate
    size(set_c,1),1);
% initialize values of c
for i = length(set_c)
    filters(i).c = set_c(i);
end
% sequence of the (indexes) tollerances selected by the controller
c_index = ones(1,steps_sim);

for t = 1:steps_sim
    %% Output of the system
    simY(:,t) = ComputeOutput(sys, simX(:,t));

    %% Controller and Filter:
    tic;
    % 1) update input of the system
    [simU(t,:), c_index(t)] = Controller(model_con, reference, simY, filters, t, con_params);

    % 2) update filters prediction of the next state
    for filt=filters
        [filt.x_pred(:,t+1), filt.V(:,:,t+1), ~, filt.P(:,:,t+1), filt.lambda(1,t+1)] = ...
            RobustKalmanFilter(model_con, ...
            filters(c_index(t)).V(:,:,t), ...      % all the filters use what the controller
            filters(c_index(t)).x_pred(:,t), ... % has determined to be the best option
            simY(:,t), simU(:,t), filt.c);
    end
    cpuT(k) = toc;

    %% Simulate the system
    simX(:,t+1) = Simulate(model_sim, simX(:,t), simU(:,t));
end
end