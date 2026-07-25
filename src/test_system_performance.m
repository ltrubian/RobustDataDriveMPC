% TEST_SYSTEM_PERFORMANCE.m
% Evaluates the overall performance of the Robust Data-Driven MPC system 
% through Monte Carlo simulations to account for parametric and additive noise.
%
% All metrics are normalized to be directly interpretable:
%   - Tracking error as % of the reference magnitude
%   - Control effort as % of the actuator range
%   - Steady-state error separated from transient behavior
%   - Constraint satisfaction as pass/fail

addpath("Controller/")
addpath("RobustKalmanFilter/")
addpath("LeastFavorableModel/")

%% Configuration
N_mc = 20;            % Number of Monte Carlo simulations
verbose = false;
model = 0;
measure_noise = true;
process_noise = true;
delta = 0.05;
offset_free = true;
steps_sim = 100;
con_params.N = 20;
con_params.L = 10;
con_params.beta = 1;
con_params.mpc = "RKF";
con_params.c_selec = "own";
set_c = [0, logspace(-6, -1, 9)];

% Reference parameters
ref_value = 5;        % Step reference magnitude
t_step = 10;          % Time step when reference jumps
t_settle = 30;        % After this step, we consider the system "settled"

if exist("osqp","class")
    con_params.options = [];
else
    con_params.options = optimoptions('quadprog', 'Display', 'off');
end

% Preallocate metric arrays
tracking_pct         = zeros(N_mc, 1);  % Tracking RMSE as % of reference
ss_error_pct         = zeros(N_mc, 1);  % Steady-state error as % of reference
control_pct          = zeros(N_mc, 1);  % Control effort as % of actuator range
estimation_pct       = zeros(N_mc, 1);  % Estimation RMSE as % of state magnitude
constraint_viol_pct  = zeros(N_mc, 1);  % Max constraint violation as % of actuator range

fprintf('Running Monte Carlo simulations (%d runs)...\n', N_mc);
tic;

for i = 1:N_mc
    rng(i); % Different seed for each run
    
    % Generate perturbed model
    [model_sim, model_con, init_con] = models(model, delta, offset_free);
    p = size(model_sim.C,1);
    n = size(model_sim.A, 1);
    m = size(model_sim.K, 2);
    
    % Reference signal
    reference = [zeros(p, t_step), ones(p, steps_sim - t_step) * ref_value];
    reference_ext = [reference, repmat(reference(:,end), 1, con_params.N)];
    
    % Actuator range
    u_range = model_con.u_max(1) - model_con.u_min(1);  % total range per input
    
    % Run simulation
    [simX, simY, trueY, simU, cpuT, filters, c_index] = LoopSimulation(...
        model_sim, model_con, steps_sim, init_con, reference_ext, set_c, measure_noise, process_noise, verbose, con_params);
    
    % --- 1. Tracking Error (% of reference) ---
    % Only measured after the step is applied (t_step onward)
    err_tracking = trueY(:, t_step+1:steps_sim) - reference(:, t_step+1:steps_sim);
    rmse_track = sqrt(mean(err_tracking(:).^2));
    tracking_pct(i) = (rmse_track / ref_value) * 100;
    
    % --- 2. Steady-State Error (% of reference) ---
    % Measured only in the settled region (after t_settle)
    err_ss = trueY(:, t_settle+1:steps_sim) - reference(:, t_settle+1:steps_sim);
    rmse_ss = sqrt(mean(err_ss(:).^2));
    ss_error_pct(i) = (rmse_ss / ref_value) * 100;
    
    % --- 3. Control Effort (% of actuator range) ---
    % RMS of each input channel, normalized by the available range
    rms_u = sqrt(mean(simU.^2, 2));   % RMS per channel [m x 1]
    control_pct(i) = mean(rms_u / u_range) * 100;
    
    % --- 4. Constraint Violation (% of actuator range) ---
    viol_min = max(0, model_con.u_min' - simU);
    viol_max = max(0, simU - model_con.u_max');
    max_viol = max([viol_min(:); viol_max(:)]);
    constraint_viol_pct(i) = (max_viol / u_range) * 100;
    
    % --- 5. Estimation Error (% of state magnitude) ---
    x_hat = zeros(n, steps_sim);
    for t = 1:steps_sim
        idx = c_index(t+1);
        x_hat(:, t) = filters(idx).x_pred(1:n, t+1);
    end
    err_est = simX(:, 2:steps_sim+1) - x_hat;
    % Normalize by the RMS of the true state (avoids division by zero at start)
    x_rms = sqrt(mean(simX(:, t_settle+1:steps_sim+1).^2, 'all'));
    rmse_est = sqrt(mean(err_est(:, t_settle:end).^2, 'all'));
    estimation_pct(i) = (rmse_est / x_rms) * 100;
    
    if mod(i, 5) == 0
        fprintf('  Completed %d / %d\n', i, N_mc);
    end
end
sim_time = toc;

%% Results Summary
fprintf('\n============================================\n');
fprintf('  PERFORMANCE REPORT  (%d Monte Carlo runs)\n', N_mc);
fprintf('  Reference = %.1f | Actuator range = [%.1f, %.1f]\n', ...
    ref_value, model_con.u_min(1), model_con.u_max(1));
fprintf('  Uncertainty level delta = %.2f\n', delta);
fprintf('============================================\n\n');

fprintf('1. TRACKING ERROR (after step, %% of reference = %.1f)\n', ref_value);
fprintf('   Mean: %5.1f%%  |  Worst-case: %5.1f%%\n', ...
    mean(tracking_pct), max(tracking_pct));
fprintf('   -> "On average, the output deviates %.1f%% from the target"\n\n', ...
    mean(tracking_pct));

fprintf('2. STEADY-STATE ERROR (after settling, t > %d, %% of reference)\n', t_settle);
fprintf('   Mean: %5.1f%%  |  Worst-case: %5.1f%%\n', ...
    mean(ss_error_pct), max(ss_error_pct));
fprintf('   -> "Once settled, the output is within %.1f%% of the target"\n\n', ...
    mean(ss_error_pct));

fprintf('3. CONTROL EFFORT (%% of actuator range = %.1f)\n', u_range);
fprintf('   Mean: %5.1f%%  |  Worst-case: %5.1f%%\n', ...
    mean(control_pct), max(control_pct));
fprintf('   -> "The controller uses %.1f%% of the available actuator capacity"\n\n', ...
    mean(control_pct));

fprintf('4. MAXIMUM CONSTRAINT VIOLATION (%% of actuator range)\n');
fprintf('   Mean: %5.1f%%  |  Worst-case: %5.1f%%\n', ...
    mean(constraint_viol_pct), max(constraint_viol_pct));
fprintf('   -> "On average, the worst constraint violation is %.1f%% of actuator range"\n\n', ...
    mean(constraint_viol_pct));

fprintf('5. STATE ESTIMATION ERROR (after settling, %% of state magnitude)\n');
fprintf('   Mean: %5.1f%%  |  Worst-case: %5.1f%%\n', ...
    mean(estimation_pct), max(estimation_pct));
fprintf('   -> "The filter estimates the true state within %.1f%% accuracy"\n\n', ...
    mean(estimation_pct));

fprintf('============================================\n');
fprintf('  Completed in %.2f seconds\n', sim_time);
fprintf('============================================\n');

%% Visualizations
figure('Name', 'Monte Carlo Performance Metrics', 'Position', [100 100 900 400]);
t = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, sprintf('Performance over %d Monte Carlo runs (\\delta = %.2f)', N_mc, delta));

nexttile;
boxchart([tracking_pct, ss_error_pct]);
xticklabels({'Overall', 'Steady-State'});
ylabel('% of reference');
title('Tracking Error');
yline(10, '--r', '10% threshold');

nexttile;
boxchart(control_pct);
xticklabels({'Effort'});
ylabel('% of actuator range');
title('Control Effort');
yline(50, '--r', '50% threshold');

nexttile;
boxchart(estimation_pct);
xticklabels({'Filter'});
ylabel('% of state magnitude');
title('Estimation Error');
yline(20, '--r', '20% threshold');

nexttile;
boxchart(constraint_viol_pct);
xticklabels({'Violation'});
ylabel('% of actuator range');
title('Max Constraint Viol.');
yline(0, '--r', '0% target');

%% Helper
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
