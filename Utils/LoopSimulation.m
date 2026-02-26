function [simX, simU, cpuT] = LoopSimulation(model_sim, model_con, ...
    steps_sim, freq_con, ini_con, reference, set_c)
%LOOPSIMULATION Simulate cloosed-loop system
%
%       <usage here>
% INPUT:
%   model_sim:  model to simulate
%   model_con:  nominal startgin model used by MPC
%   steps_sim:  number of step to simulate
%   freq_con:   how many steps the controller updates the input (MPC
%               frequency
%   init_con:   initial condition
%   reference:  reference signal
%   set_c:      set of hyperparamter 'c' to choose from
%
% OUTPUT:
%   simX: simulated states
%   simU: controlled input
%   cpuT: cpu time of the controller
%
% DEV-STATUS:

arguments
    model_sim   (1,1) struct
    model_con   (1,1) struct
    steps_sim   (1,1) double {mustBeInteger(steps_sim)}
    freq_con    (1,1) double {mustBeInteger(freq_con)}
    ini_con     (:,:) double
    reference   (:,:) double
    set_c       (1,:) double
end

n = size(model_sim.A,1);
p = size(model_sim.C,1);
m = size(model_sim.B,2);

% state and control trajectories
simX = zeros(steps_sim + 1, n);
simY = zeros(steps_sim + 1, p);
simU = zeros(steps_sim / freq_con, m);
cpuT = zeros(size(simU,1), 1);
simX(1,:) = ini_con;

% controller index
k = 1;

for i = 1:steps_sim
    if mod(i, freq_con) == 0
        tic;
        % update controller
        % simU(k,:) = Controller(model_con, simX(i-1,:), reference(k,:), set_c)
        cpuT(k) = toc;
        k = k + 1;
    end
    % simulate the system
    % [simX(i,:), simY(i,:)] = simulation(model_sim, simX(i-1,:), simU(k-1,:));
end
end