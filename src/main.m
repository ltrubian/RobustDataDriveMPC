addpath("utils/")

%% ----------------------- CONFIG -----------------------

% Enable terminal output during execution
verbose = false;

% True model simulation noises config
noise_config = struct( ...
    "measure_noise", true, ...  % if false set measure noise to zero
    "process_noise", true, ...  % if false set process noise to zero
    "delta_process", 0.05, ...  % Noise matrices diagonal elements
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

try
    % brief check on eventual saturation problem of actuator at steady-state
    uss  = [model_sim.A - eye(n), model_sim.K; model_sim.C, zeros(p,m)] \ [zeros(n,1); reference(:,end)];
    uss = uss(n+1:end);
    % stochastic analysis
    % static input gain
    F = dlqr(model_sim.A, model_sim.K, [model_nom.weights.Q, zeros(p,2); zeros(2,p) eye(2)], model_nom.weights.R);
    % static kalman gain
    L = dlqe(model_sim.A, eye(n), model_sim.C, model_sim.B*model_sim.B', model_sim.D*model_sim.D', model_sim.B*model_sim.D');
    % closed-loop matrices
    A = [model_sim.A - model_sim.K * F, model_sim.K * F; zeros(n,n), model_sim.A - model_sim.A* L * model_sim.C ];
    B = [model_sim.B; model_sim.B - L * model_sim.D];
    S1 = dlyap(A, B*B');
    S = [-F, F] * S1 * [-F, F]';
    % stochastic bounds
    [uss + 3 * sqrt(diag(S))]'
    [uss - 3 * sqrt(diag(S))]'
catch
end
% model_nom.u_max = [uss + 10 * sqrt(diag(S))]';
% model_nom.u_min = [uss - 10 * sqrt(diag(S))]';

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