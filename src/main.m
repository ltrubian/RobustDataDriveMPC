addpath("utils/")

%% ----------------------- CONFIG -----------------------

% Enable terminal output during execution
verbose = false;

% True and nominal models noises config
noise_config = struct( ...
    "measure_noise", true, ...  % enable measure noise in true model simulation
    "process_noise", true, ...  % enable process noise in true model simulation
    "delta_process", 0.05, ...  % nominal model noise matrices diagonal elements
    "delta_measure", 0.05  ...  %   "
);

% MPC Controller configuration
MPC_config = struct( ...
    "offset_free", true, ...
    "u_min",       0, ...
    "u_max",       5, ...
    "x_min",       0, ...
    "x_max",       30, ...
    "Q",           1, ...
    "Pf",          1, ...
    "R",           1 ...
);
% Note
% offset_free: introduces fictitious constant disturbances in the nominal
% model allowing the MPC to compensate for modelling errors

% Algorithm-side configuration
con_params = struct( ...
    "N",        20, ...      % prediction horizon of MPC
    "L",        10, ...      % time windows toward the past for estimation
    "beta",     1, ...       % forgetting factor
    "mpc",      "RKF", ...   % MPC strategy ("RKF" or "LFM")
    "c_selec",  "own", ...   % filter selection mode ("best" or "own")
    "options",  optimoptions('quadprog', ...
                    'OptimalityTolerance', 1e-6, ...
                    'StepTolerance', 1e-6, ...
                    'ConstraintTolerance', 1e-6, ...
                    'Display', 'off') ...
);

if exist("osqp","class")
    con_params.options = [];
else
    warning("consider installing osqp solver for faster execution")
end

% steps_sim: number of step to simulate
steps_sim = 200;

% set_c: set of hyperparamter 'c' to choose from
set_c = [0, logspace(-6, -1, 9)];

% ----------------------- END-CONFIG -----------------------

[model_sim, model_nom, init_con] = model_config(noise_config, MPC_config);

n = size(model_sim.A,1);        % state real world
m = size(model_sim.K,2);        % input real world
p = size(model_sim.C,1);        % output real world
r = size(model_nom.A,1) - n;    % fictitious disturbances (if introduced)

% reference: reference signal
reference = [ones(p, 10) * 10, ones(p, steps_sim-10) * 15];

% update reference: last value is repeated so that the controller has enough preview
reference = [reference, repmat(reference(:,end), 1, con_params.N)];

%% SIMULATION LOOP
[simX, simY, trueY, simU, cpuT, filters, c_index] = closed_loop_simulation(model_sim, model_nom, ...
    steps_sim, init_con, reference, set_c, noise_config.measure_noise, noise_config.process_noise, verbose, ...
    con_params);

%% RECAP simulation params and controller strategy
fprintf("=== simulaition params  ===\n")
fprintf("process noise: %s\nmeasure noise: %s\ndelta process: %.2f\ndelta measure: %.2f\n", ...
    string(noise_config.process_noise), string(noise_config.measure_noise), ...
    noise_config.delta_process, noise_config.delta_measure)
fprintf("=== controller strategy ===\n")
fprintf("mpc strategy: %s\n c selection: %s\n offset-free: %s\n", ...
    con_params.mpc, con_params.c_selec, string(MPC_config.offset_free))
%% Plot and analysis
plotting;