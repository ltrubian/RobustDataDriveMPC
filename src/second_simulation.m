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

[model_sim, model_con, init_con] = models(2, delta*(1-debug));
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
    x_hat(:,t) = filters(c_index(t)).x_pred(1:n,t);
    y_hat(:,t) = model_con.C * x_hat(:,t);
end

% Plot results
figure;

ax1=subplot(3,1,1);
plot((simX(1,1:end-1)-x_hat(1,:)).^2);
hold on;
plot((simX(2,1:end-1)-x_hat(2,:)).^2);
legend('er_1', 'er_2');
title('True VS Estimated States');

ax2=subplot(3,1,2);
plot(simY(1,:));
hold on;
plot(y_hat(1,:));
hold on;
plot(reference(1,1:steps_sim));
legend('True Output','RKF Estimate');
title('True VS Estimated Output');

ax3=subplot(3,1,3);
stairs(simU');
title('Ingresso di Controllo u');

xlim([ax1,ax2,ax3],[1, steps_sim])