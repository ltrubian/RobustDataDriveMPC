addpath("Controller/")
addpath("RobustKalmanFilter/")
rng(1)
%% DEFINITION OF VARIABLES FOR THE SIMULATION
%   model_sim:  model to simulate
model_sim.A = [1.1 1; 0 1];
model_sim.B = [0.5; 1];
model_sim.C = [1 0];
model_sim.D = [0.1];

n = size(model_sim.A,1);
p = size(model_sim.C,1);
m = size(model_sim.B,2);

%   model_con:  nominal startgin model used by MPC
model_con.A = [1.15 1; 0.01 1];
model_con.B = [0.48; 1.05];
model_con.C = [0.99 0];
model_con.D = 0;
model_con.u_min = -2; model_con.u_max = 2;
model_con.x_min = [-inf; -inf];
model_con.x_max = [+inf; +inf];
model_con.weights.Q = 1;
model_con.weights.Pf = 1;
model_con.weights.R = 0.01;

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
con_params.N = 10;
reference = [reference,repmat(reference(end),1,con_params.N)];
%       L:      time windows relevant for estimation
con_params.L = 10;
%       beta:   forgetting factor
con_params.beta = 0.9;

%% SIMULATION OF THE WHOLE SYSTEM
[simX, simY, simU, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
    steps_sim, init_con, reference, set_c', L=con_params.L, N=con_params.N, beta=con_params.beta);

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
% (and the filters)
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

