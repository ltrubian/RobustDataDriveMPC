addpath("Controller/")
addpath("RobustKalmanFilter/")
addpath("LeastFavorableModel/")
rng(1)
verbose = true;

%% DEFINITION OF VARIABLES FOR THE SIMULATION
%   model_sim:  true model to simulate
model_sim.A = [1.1 1; 0 1];         % state -> state
model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];   % noise -> state
model_sim.C = [1 0];                % state -> output
model_sim.D = [0.1, 0.05, 0.01];          % noise -> output

model_sim.K = [0.5; 1];             % input -> state
% model_sim.J = [0.1; 0.05];          % input -> output

n = size(model_sim.A,1);
p = size(model_sim.C,1);
m = size(model_sim.K,2);

% Struct containing all the gains for noises/disturbances
delta = 0.01;      % model perturbation gain

% "DEBUG MODE": if True set all the noise/perturbation gains to 0
% Use to check if the MPC controller works in ideal conditions
debug = false;

if debug
    delta = 0;
end

%   model_con:  nominal (perturbed) model used by MPC controller. The
%   perturbation of each entry is the product of the gain delta and a
%   random matrix with compatible sie
model_con.A = model_sim.A + delta * randn(size(model_sim.A));
model_con.B = model_sim.B + delta * randn(size(model_sim.B));
model_con.C = model_sim.C + delta * randn(size(model_sim.C));
model_con.D = model_sim.D + delta * randn(size(model_sim.D));

model_con.K = model_sim.K + delta * randn(size(model_sim.K));
% model_con.J = model_sim.J + delta * randn(size(model_sim.J));

% MPC config
model_con.u_min = -2;
model_con.u_max = 2;
model_con.x_min = [-inf; -inf];
model_con.x_max = [+inf; +inf];
model_con.weights.Q = 1;
model_con.weights.Pf = 1;
model_con.weights.R = 0.1;

%   steps_sim:  number of step to simulate
steps_sim = 50;

%   init_con:   initial condition
init_con = [1; 0];

%   reference:  reference signal
reference = ones(1, steps_sim) * 5;
% time = 1:steps_sim;
% reference = sin(0.1*time);

%   set_c:      set of hyperparamter 'c' to choose from
set_c = [0, logspace(-6, -1, 9)];

% NAMED-VALUE INPUTS:
%   con_params:
%       N:      prediction horizon of MPC
con_params.N = 20;
% update reference: last value is repeated so that the controller has
% always enough preview
reference = [reference,repmat(reference(end),1,con_params.N)];
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
con_params.mpc = "RKF-ext";
% con_params.mpc = "RKF";
% con_params.mpc = "LFM";

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
[simX, simY, simU, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
    steps_sim, init_con, reference, set_c, debug, verbose, ...
    con_params);

%% report and analysis
% OUTPUT OF THE SIMULATION:
%   simX:       simulated states
%   simY:       simulated output
%   simU:       controlled input
%   cpuT:       cpu time of the controller
%   filters:    struct with the dynamincs of the set of filters
%   c_index:    the sequence of c's selected by the controller
%
% COMPUTED/EXTRACTED VALUES:
%   x_hat:      state estimation made by the (combination of) controller and
%               filters and used by the first to produce the input and the
%               latter estimate the next state

% the struct filters holds all the prediction made by all the filters and
% we need to extract the ones that were actually used by the controller
x_hat = zeros(n, steps_sim);
y_hat = zeros(p, steps_sim);
for t=1:steps_sim
    x_hat(:,t) = filters(c_index(t)).x_pred(:,t);
    y_hat(:,t) = model_con.C * x_hat(:,t);
end

% Plot results
figure;

ax1=subplot(3,1,1);
plot((simX(1,2:end)-x_hat(1,:)).^2);
hold on;
plot((simX(2,2:end)-x_hat(2,:)).^2);
legend('er_1', 'er_2');
title('True VS Estimated States');

ax2=subplot(3,1,2);
plot(simY(1,:));
hold on;
plot(y_hat(1,:));
hold on;
plot(reference(1:steps_sim));
legend('True Output','RKF Estimate');
title('True VS Estimated Output');

ax3=subplot(3,1,3);
stairs(simU(1,:));
title('Ingresso di Controllo u');

xlim([ax1,ax2,ax3],[1, steps_sim])
