addpath("../Controller/")
addpath("../RobustKalmanFilter/")
addpath("../LeastFavorableModel/")
addpath("../")
rng(1)
verbose = false;
% the same true model and nominal model (used by the controller) is used
% for all the simulations. Each simulation varies from the others on the
% random noise.
% the following parameters are crucial to decide how long this script will
% take:
%   - steps_sim     : how many steps each simulation should last
%   - n_simul       : how many different simulation will be run for each
%                     controller type

%% DEFINITION OF VARIABLES FOR THE SIMULATION
%   model_sim:  true model to simulate
model_sim.A = [1.1 1; 0 1];         % state -> state
model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];   % noise -> state
model_sim.C = [1 0];                % state -> output
model_sim.D = [0.1, 0.05, 0.01];          % noise -> output

model_sim.K = [0.5 0; 1 0.1];             % input -> state
% model_sim.J = [0.1; 0.05];          % input -> output

n = size(model_sim.A,1);
p = size(model_sim.C,1);
m = size(model_sim.K,2);

% Struct containing all the gains for noises/disturbances
delta = 0.05;      % model perturbation gain

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
model_con.u_min = [-2, -2];
model_con.u_max = [2, 2];
model_con.x_min = [-inf; -inf];
model_con.x_max = [+inf; +inf];
model_con.weights.Q = 1;
model_con.weights.Pf = 1;
model_con.weights.R = 0.1*eye(2);

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
con_params.beta = 0.95;
%       mpc:    which strategy to use the MPC
%               RKF-ext: exted the model to the N time horizon and make RKF
%                        to that extended model (just starting point x0 is
%                        given to MPC)
%               RKF:     compute the RKF on the nominal model (just
%                        starting point x0 is given to MPC)
%               LFM:     the time-varying LFM is computed and used for the
%                        prediction x0 (LFM model and x0 are given to MPC)
%
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

%% MULTIPLE SIMULATIONS OF THE SAME SYSTEM
% number of simulations to run: for each simulation one unique seed is used
% for both type of controllers. Increasing this value increases the
% execution time of this script and the accuracy of the results
n_simul = 10;
% starting seed: all the simulations are done with the seed <i + s_simul>.
% in order to make different runs of the script you need to vary this one
s_simul = 5000000;

% struct to collect errors along time for Least-Favorable Model
lfm.ex_pred = NaN(n, steps_sim, n_simul);
lfm.ey_pred = NaN(p, steps_sim, n_simul);
lfm.ey      = NaN(p, steps_sim, n_simul);
% and for the Robust Kalman Filter approach on the extended state
% [x_k, x_k+1, ..., x_k+N]
rkf.ex_pred = NaN(n, steps_sim, n_simul);
rkf.ey_pred = NaN(p, steps_sim, n_simul);
rkf.ey      = NaN(p, steps_sim, n_simul);

% total time of execution of these long computations
total_time = 0;
for i=1:n_simul
    % if i>2; break; end
    loops = tic;
    % test the LFM formulation of the controller
    con_params.mpc = "LFM";
    % set the random seed and simulate
    rng(i+s_simul);
    [simX, simY, ~, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
        steps_sim, init_con, reference, set_c, debug, verbose, ...
        con_params);
    % collect preditions of the filters
    x_hat = zeros(n, steps_sim);
    y_hat = zeros(p, steps_sim);
    for t=1:steps_sim
        x_hat(:,t) = filters(c_index(t)).x_pred(1:n,t);
        y_hat(:,t) = model_con.C * x_hat(:,t);
    end
    % compute and store errors
    lfm.ex_pred(:,:,i) = (simX(:,1:end-1) - x_hat).^2;
    lfm.ey_pred(:,:,i) = (simY(:,1:end-1) - y_hat).^2;
    lfm.ey(:,:,i)      = (simY(:,1:end-1) - reference(1:steps_sim)).^2;

    fprintf("sim: %2d/%2d, type: %7s, time: %.2e err_x1: %2.5f\n",...
        i, n_simul, con_params.mpc, sum(cpuT), mean(lfm.ex_pred(1,:,i)))

    % test the LFM formulation of the controller
    con_params.mpc = "RKF";
    % set the random seed and simulate
    rng(i+s_simul);
    [simX, simY, ~, cpuT, filters, c_index] = LoopSimulation(model_sim, model_con, ...
        steps_sim, init_con, reference, set_c, debug, verbose, ...
        con_params);
    % collect preditions of the filters
    x_hat = zeros(n, steps_sim);
    y_hat = zeros(p, steps_sim);
    for t=1:steps_sim
        x_hat(:,t) = filters(c_index(t)).x_pred(1:n,t);
        y_hat(:,t) = model_con.C * x_hat(:,t);
    end
    % compute and store errors
    rkf.ex_pred(:,:,i) = (simX(:,1:end-1) - x_hat).^2;
    rkf.ey_pred(:,:,i) = (simY(:,1:end-1) - y_hat).^2;
    rkf.ey(:,:,i)      = (simY(:,1:end-1) - reference(1:steps_sim)).^2;

    fprintf("sim: %2d/%2d, type: %7s, time: %.2e err_x1: %2.5f\n",...
        i, n_simul, con_params.mpc, sum(cpuT), mean(rkf.ex_pred(1,:,i)))

    % estimate of the remaining time of execution
    current_time = toc(loops);
    total_time = current_time + total_time;
    fprintf("  current: %.2f \t remaining: %.2f \t total: %.2f\n", ...
        current_time, total_time / i * (n_simul - i), total_time)
end
%% plots
% considering the error at time t in the simulation k <e_t^k>
% the plot represents the mean of the error at time t for all the
% simulation so the graph will show for each t <mean_k(e_t^k)>
figure(Name="Errors evolution")
stadard_dev = true;
ti = 30:steps_sim;
H = rgb2hex(orderedcolors("gem"));

% Prediction errors on state 1
ax0=subplot(2,2,1);
% compute and plot mean across multiple runs 
mean_lfm = mean(lfm.ex_pred(1,ti,:),3, "omitnan");
mean_rkf = mean(rkf.ex_pred(1,ti,:),3, "omitnan");
plot(ti, mean_lfm, Color=H(1)); hold on;
plot(ti, mean_rkf, Color=H(2)); hold on;
% compute and plot standard deviations
if stadard_dev
    stde_lfm = std(lfm.ex_pred(1,ti,:),0,3, "omitnan");
    stde_rkf = std(rkf.ex_pred(1,ti,:),0,3, "omitnan");
    plot(ti, mean_lfm+stde_lfm,ti, mean_lfm-stde_lfm, Color=H(1), LineStyle="--" ); hold on;
    plot(ti, mean_rkf+stde_rkf,ti, mean_rkf-stde_rkf, Color=H(2), LineStyle="--" ); hold off;
end
legend('lfm', 'rkf');
title('errors on 1st state');

% Prediction errors on state 1
ax1=subplot(2,2,2);
% compute and plot mean across multiple runs
mean_lfm = mean(lfm.ex_pred(2,ti,:),3, "omitnan");
mean_rkf = mean(rkf.ex_pred(2,ti,:),3, "omitnan");
plot(ti, mean_lfm, Color=H(1)); hold on
plot(ti, mean_rkf, Color=H(2)); hold on
% compute and plot standard deviations
if stadard_dev
    stde_lfm = std(lfm.ex_pred(2,ti,:),0,3, "omitnan");
    stde_rkf = std(rkf.ex_pred(2,ti,:),0,3, "omitnan");
    plot(ti, mean_lfm+stde_lfm,ti, mean_lfm-stde_lfm, Color=H(1), LineStyle="--" ); hold on;
    plot(ti, mean_rkf+stde_rkf,ti, mean_rkf-stde_rkf, Color=H(2), LineStyle="--" ); hold off;
end
legend('lfm', 'rkf');
title('errors on 2nd state');

% Prediction errors on output 1
ax2=subplot(2,2,3);
% compute and plot mean across multiple runs
mean_lfm = mean(lfm.ey_pred(1,ti,:),3, "omitnan");
mean_rkf = mean(rkf.ey_pred(1,ti,:),3, "omitnan");
plot(ti, mean_lfm, Color=H(1)); hold on
plot(ti, mean_rkf, Color=H(2)); hold on
% compute and plot standard deviations
if stadard_dev
    stde_lfm = std(lfm.ey_pred(1,ti,:),0,3, "omitnan");
    stde_rkf = std(rkf.ey_pred(1,ti,:),0,3, "omitnan");
    plot(ti, mean_lfm+stde_lfm,ti, mean_lfm-stde_lfm, Color=H(1), LineStyle="--" ); hold on;
    plot(ti, mean_rkf+stde_rkf,ti, mean_rkf-stde_rkf, Color=H(2), LineStyle="--" ); hold off;
end
legend('lfm', 'rkf');
title('errors on output prediction');

% Tracking error
ax3=subplot(2,2,4);
% compute and plot mean across multiple runs
mean_lfm = mean(lfm.ey(1,ti,:),3, "omitnan");
mean_rkf = mean(rkf.ey(1,ti,:),3, "omitnan");
plot(ti, mean_lfm, Color=H(1)); hold on
plot(ti, mean_rkf, Color=H(2)); hold on
% compute and plot standard deviations
if stadard_dev
    stde_lfm = std(lfm.ey(1,ti,:),0,3, "omitnan");
    stde_rkf = std(rkf.ey(1,ti,:),0,3, "omitnan");
    plot(ti, mean_lfm+stde_lfm,ti, mean_lfm-stde_lfm, Color=H(1), LineStyle="--" ); hold on;
    plot(ti, mean_rkf+stde_rkf,ti, mean_rkf-stde_rkf, Color=H(2), LineStyle="--" ); hold off;
end
legend('lfm', 'rkf');
title('tracking error Y');

xlim([ax0,ax1,ax2,ax3],[ti(1), ti(end)])

%% comprehensive results
% the means are first evaluated along the simulations (the same as before)
% and then are evaluated along the time. In theory this should not change
% anything, but in practice this is different when some simulation
% interrupt before the end.
ti = 30:steps_sim;
tab = table(zeros(4,1), zeros(4,1), zeros(4,1), zeros(4,1), ...
    'RowNames',["ex1_pred","ex2_pred","ey_pred","ey"], ...
    'VariableNames', ["LFM", "RKF", "LFM (std)", "RKF (std)"]);

tab("ex1_pred","LFM") = {mean(mean(lfm.ex_pred(1,ti,:),2),3, "omitnan")};
tab("ex2_pred","LFM") = {mean(mean(lfm.ex_pred(2,ti,:),2),3, "omitnan")};
tab("ey_pred","LFM")  = {mean(mean(lfm.ey_pred(:,ti,:),2),3, "omitnan")};
tab("ey","LFM")       = {mean(mean(lfm.ey(:,ti,:),2),3, "omitnan")};

tab("ex1_pred","RKF") = {mean(mean(rkf.ex_pred(1,ti,:),2),3, "omitnan")};
tab("ex2_pred","RKF") = {mean(mean(rkf.ex_pred(2,ti,:),2),3, "omitnan")};
tab("ey_pred","RKF")  = {mean(mean(rkf.ey_pred(:,ti,:),2),3, "omitnan")};
tab("ey","RKF")       = {mean(mean(rkf.ey(:,ti,:),2),3, "omitnan")};

tab("ex1_pred","LFM (std)") = {mean(std(lfm.ex_pred(1,ti,:),0,2),3, "omitnan")};
tab("ex2_pred","LFM (std)") = {mean(std(lfm.ex_pred(2,ti,:),0,2),3, "omitnan")};
tab("ey_pred","LFM (std)")  = {mean(std(lfm.ey_pred(:,ti,:),0,2),3, "omitnan")};
tab("ey","LFM (std)")       = {mean(std(lfm.ey(:,ti,:),0,2),3, "omitnan")};

tab("ex1_pred","RKF (std)") = {mean(std(rkf.ex_pred(1,ti,:),0,2),3, "omitnan")};
tab("ex2_pred","RKF (std)") = {mean(std(rkf.ex_pred(2,ti,:),0,2),3, "omitnan")};
tab("ey_pred","RKF (std)")  = {mean(std(rkf.ey_pred(:,ti,:),0,2),3, "omitnan")};
tab("ey","RKF (std)")       = {mean(std(rkf.ey(:,ti,:),0,2),3, "omitnan")};

disp(tab)