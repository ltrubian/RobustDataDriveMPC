addpath("utils/")
rng(1)
verbose = true;

%% DEFINITION OF VARIABLES FOR THE SIMULATION
% if false set measure noise to zero
measure_noise = true;
% if false set process noise to zero
process_noise = true;
% Noise matrices diagonal elements
delta_process = 0.05;
delta_measure = 0.05;
% offset_free: introduces fictitious constant disturbances in the nominal
% model allowing the MPC to compensate for offset
offset_free = true;

[model_sim, model_nom, init_con] = models(9, delta_process, delta_measure, offset_free);

n = size(model_sim.A,1);        % state real world
m = size(model_sim.K,2);        % input real world
p = size(model_sim.C,1);        % output real world
r = size(model_nom.A,1) - n;    % fictitious disturbances (if introduced)

%   steps_sim:  number of step to simulate
steps_sim = 200;

%   reference:  reference signal
reference = [ones(p, 10) * 10, ones(p, steps_sim-10) * 15];
try
    % % brief check on eventual saturation problem of actuator at steady-state
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
%               RKF:     compute the RKF on the nominal model (just
%                        starting point x0 is given to MPC)
%               LFM:     the time-varying LFM is computed and used for the
%                        prediction x0 (LFM model and x0 are given to MPC)
choice = 1;
switch choice
    case 1
        con_params.mpc = "RKF";
    case 2
        con_params.mpc = "LFM";
end
%       c_selection: which V and x_pred the filters will use.
%                    each filter uses
%                    best: the best x_pred (and V) of previous round
%                    own:  its own x_pred and V
%       WARNING: in the RKF approach the "best" selection is useless since
%                the dynamics of the whole algorithm will always goes for
%                the first c of the list
choice = 2;
switch choice
    case 1
        con_params.c_selec = "best";
    case 2
        con_params.c_selec = "own";
end

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
[simX, simY, trueY, simU, cpuT, filters, c_index] = closed_loop_simulation(model_sim, model_nom, ...
    steps_sim, init_con, reference, set_c, measure_noise, process_noise, verbose, ...
    con_params);

%% RECAP simulation params and controller strategy
fprintf("=== simulaition params  ===\n")
fprintf("process noise: %s\nmeasure noise: %s\n\tdelta_process: %.2f\n\tdelta_measure: %.2f\n", ...
    string(process_noise), string(measure_noise), delta_process, delta_measure)
fprintf("=== controller strategy ===\n")
fprintf("mpc strategy: %s\n c selection: %s\n offset-free: %s\n", ...
    con_params.mpc, con_params.c_selec, string(offset_free))
%% Plot and analysis
plotting;