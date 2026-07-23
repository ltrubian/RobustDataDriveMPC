addpath("Controller/")
addpath("RobustKalmanFilter/")
addpath("LeastFavorableModel/")
rng(1)
verbose = true;

%% DEFINITION OF VARIABLES FOR THE SIMULATION
% "DEBUG MODE": if True set all the noise/perturbation gains to 0
% Use to check if the MPC controller works in ideal conditions
debug = false;
% Struct containing all the gains for noises/disturbances
delta = 0.05;      % model perturbation gain
% offset_free
offset_free = true;

[model_sim, model_con, init_con] = models(2, delta, offset_free);

n = size(model_sim.A,1);
m = size(model_sim.K,2);
p = size(model_sim.C,1);

%   steps_sim:  number of step to simulate
steps_sim = 100;

%   reference:  reference signal
reference = [zeros(p, 10), ones(p, steps_sim-10) * 5];

%   set_c:      set of hyperparamter 'c' to choose from
set_c = [0, logspace(-6, -1, 9)];

% NAMED-VALUE INPUTS:
%   con_params:
%       N:      prediction horizon of MPC
con_params.N = 20;
% update reference: last value is repeated so that the controller has
% always enough preview
reference = [reference,repmat(reference(:,end),1,con_params.N)];
%       L:      time windows toward the past for estimation
con_params.L = 10;
%       beta:   forgetting factor
con_params.beta = 1;
%       mpc:    which strategy to use the MPC
%               RKF-ext: exted the model to the N time horizon and make RKF
%                        to that extended model (just starting point x0 is
%                        given to MPC)
%               RKF:     compute the RKF on the nominal model (just
%                        starting point x0 is given to MPC)
%               LFM:     the time-varying LFM is computed and used for the
%                        prediction x0 (LFM model and x0 are given to MPC)
% con_params.mpc = "RKF-ext";
con_params.mpc = "RKF";
con_params.mpc = "LFM";

con_params.options = optimoptions('quadprog', ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6, ...
    'Display', 'off');

if exist("osqp","class")
    con_params.options = [];
else
    warning("consider installing osqp solver for faster execution")
end

%% SIMULATION OF THE WHOLE SYSTEM
[simX, simY, trueY, simU, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
    steps_sim, init_con, reference, set_c, debug, verbose, ...
    con_params);

%% Plot and analysis
Plotting