addpath("Controller/")
addpath("RobustKalmanFilter/")
rng(1)
%% DEFINITION OF VARIABLES FOR THE SIMULATION
%   model_sim:  true model to simulate
model_sim.A = [1.1 1; 0 1];
model_sim.B = [0.5; 1];
model_sim.C = [1 0];
model_sim.D = [0.1];

n = size(model_sim.A,1);
p = size(model_sim.C,1);
m = size(model_sim.B,2);

%   delta: perturbation gain. Used to tune the magnitude of the
%   perturbation
delta = 0.15;

%   sign: matrix of 1 and -1 to make the perturbation sign random
A_sign = 2 * randi([0, 1], n) - 1;
B_sign = 2 * randi([0, 1], n, m) - 1;
C_sign = 2 * randi([0, 1], p, n) - 1;
D_sign = 2 * randi([0, 1], p, m) - 1;

%   model_con:  nominal (perturbed) model used by MPC controller. The
%   perturbation of each entry is the product of the gain delta, 1 or -1 to
%   make the error casual and the original entry
model_con.A = model_sim.A + delta*A_sign.*model_sim.A;
model_con.B = model_sim.B + delta*B_sign.*model_sim.B;
model_con.C = model_sim.C + delta*C_sign.*model_sim.C;
model_con.D = model_sim.D + delta*D_sign.*model_sim.D;

model_con.u_min = -2; model_con.u_max = 2;
model_con.x_min = [-inf; -inf];
model_con.x_max = [+inf; +inf];
model_con.weights.Q = 1;
model_con.weights.Pf = 1;
model_con.weights.R = 0.1;

%   steps_sim:  number of step to simulate
steps_sim = 300;

%   init_con:   initial condition
init_con = [1; 0];

%   reference:  reference signal
reference = ones(1, steps_sim) * 5;

%   set_c:      set of hyperparamter 'c' to choose from
set_c = logspace(-6, -3, 10);

% NAMED-VALUE INPUTS:
%   con_params:
%       N:      prediction horizon of MPC
con_params.N = 20;
% update reference: last value is repeated so that the controller has
% always enough preview
reference = [reference,repmat(reference(end),1,con_params.N)];
%       L:      time windows relevant for estimation
con_params.L = 10;
%       beta:   forgetting factor
con_params.beta = 0.9;
%       lfm:    apply Least-Favorable Model (true/false)
con_params.lfm = true;
%       steps:  1 to combine estimation and controller;
%               2 to have estimation minimizing before and then controller
con_params.steps = 1;
%       mpc:    which matlab function to use for the MPC controller
%               fmincon: more readble but slow (not suitable for big N)
%               quadprog: fast quadratic solver for sparse mpc
%                         implemntation
con_params.mpc = "quadprog";

%% SIMULATION OF THE WHOLE SYSTEM
[simX, simY, simU, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
    steps_sim, init_con, reference, set_c', ...
    L=con_params.L, N=con_params.N, beta=con_params.beta, ...
    lfm=con_params.lfm, steps=con_params.steps, mpc=con_params.mpc);

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

subplot(3,1,1);
plot(simX(1,:));
hold on;
plot(x_hat(1,:));
legend('True State','RKF Estimate');
title('True VS Estimated States');

subplot(3,1,2);
plot(simY(1,:));
hold on;
plot(y_hat(1,:));
hold on;
plot(reference(1:steps_sim));
legend('True Output','RKF Estimate');
title('True VS Estimated Output');

subplot(3,1,3);
stairs(simU(1,:));
title('Ingresso di Controllo u');

