% TEST_SYSTEM_PERFORMANCE.m
% Evaluates the overall performance of the Robust Data-Driven MPC system
% through Monte Carlo simulations over both parametric uncertainty (random
% perturbations of the true plant matrices A, K) and additive noise
% (process and measurement noise realizations).
%
% All metrics are normalized to be directly interpretable:
%   - Tracking error as % of the reference level
%   - Control effort as % of the actuator range
%   - Steady-state error separated from transient behavior
%   - Constraint satisfaction as pass/fail

addpath("utils/")

%% Configuration
N_mc = 20;            % Number of Monte Carlo simulations
verbose = false;

% True and nominal models noises config
noise_config = struct( ...
    "measure_noise", true, ...
    "process_noise", true, ...
    "delta_process", 0.05, ...
    "delta_measure", 0.05  ...
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

% Simulation length
steps_sim = 100;

% Parametric uncertainty: element-wise std deviation of multiplicative
% perturbation applied to the true plant matrices A and K each run.
% Set to 0 to disable parametric variation (noise-only Monte Carlo).
param_pert_scale = 0.01;

% Algorithm-side configuration
con_params = struct( ...
    "N",        20, ...
    "L",        10, ...
    "beta",     1, ...
    "mpc",      "LFM", ...
    "c_selec",  "own", ...
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

% set_c: set of hyperparameter 'c' to choose from
set_c = [0, logspace(-6, -1, 9)];

% Reference parameters
t_step = 10;          % Time step when reference jumps
t_settle = 50;        % After this step, we consider the system "settled"

% Build a temporary model to determine output dimension p for reference
[tmp_sim, ~, ~] = model_config(noise_config, MPC_config);
p = size(tmp_sim.C, 1);
clear tmp_sim

% Reference signal
reference = [ones(p, t_step) * 10, ones(p, steps_sim-t_step) * 15];
reference_ext = [reference, repmat(reference(:,end), 1, con_params.N)];

% Normalization: use the post-step reference level (not the step size)
ref_level = reference(1, end);  % = 15

% Preallocate metric arrays (_pct means percent)
tracking_pct         = zeros(N_mc, 1);  % Tracking RMSE as % of reference level
ss_error_pct         = zeros(N_mc, 1);  % Steady-state error as % of reference level
control_pct          = zeros(N_mc, 1);  % Control effort as % of actuator range
estimation_pct       = zeros(N_mc, 1);  % Estimation RMSE as % of output magnitude
constraint_viol_pct  = zeros(N_mc, 1);  % Max constraint violation as % of actuator range

fprintf('Running Monte Carlo simulations (%d runs)...\n', N_mc);
fprintf('  Parametric perturbation scale: %.1f%%\n', param_pert_scale*100);

sim_time = tic;
for i = 1:N_mc
    rng(i); % Different seed for each run

    % Generate base model from deterministic config
    [model_sim, model_nom, init_con] = model_config(noise_config, MPC_config);
    n = size(model_sim.A, 1);
    m = size(model_sim.K, 2);
    p = size(model_sim.C, 1);
    r = size(model_nom.A, 1) - n;    % fictitious disturbances (if introduced)

    % --- Parametric perturbation of the true plant ---
    % Each run sees a different realization of the true A and K matrices,
    % modelling the fact that the real plant is only approximately known.
    if param_pert_scale > 0
        max_attempts = 50;
        stable = false;
        for attempt = 1:max_attempts
            A_pert = model_sim.A .* (1 + param_pert_scale * randn(n));
            K_pert = model_sim.K .* (1 + param_pert_scale * randn(n, m));
            if max(abs(eig(A_pert))) < 1
                stable = true;
                break;
            end
        end
        if ~stable
            warning('Run %d: could not find stable perturbation in %d attempts, using unperturbed model', ...
                i, max_attempts);
        else
            model_sim.A = A_pert;
            model_sim.K = K_pert;
        end
    end

    % Actuator range
    u_range = model_nom.u_max(1) - model_nom.u_min(1);

    % Run simulation
    [simX, simY, trueY, simU, cpuT, RKFs, c_best] = closed_loop_simulation(...
        model_sim, model_nom, steps_sim, init_con, reference_ext, set_c, ...
        noise_config.measure_noise, noise_config.process_noise, verbose, con_params);

    % --- 1. Tracking Error (% of reference level) ---
    % Only measured after the step is applied (t_step onward)
    err_tracking = trueY(:, t_step+1:steps_sim) - reference(:, t_step+1:steps_sim);
    rmse_track = sqrt(mean(err_tracking(:).^2));
    tracking_pct(i) = (rmse_track / ref_level) * 100;

    % --- 2. Steady-State Error (% of reference level) ---
    % Measured only in the settled region (after t_settle)
    err_ss = trueY(:, t_settle+1:steps_sim) - reference(:, t_settle+1:steps_sim);
    rmse_ss = sqrt(mean(err_ss(:).^2));
    ss_error_pct(i) = (rmse_ss / ref_level) * 100;

    % --- 3. Control Effort (% of actuator range) ---
    rms_u = sqrt(mean(simU.^2, 2));
    control_pct(i) = mean(rms_u / u_range) * 100;

    % --- 4. Constraint Violation (% of actuator range) ---
    viol_min = max(0, model_nom.u_min - simU);
    viol_max = max(0, simU - model_nom.u_max);
    max_viol = max([viol_min(:); viol_max(:)]);
    constraint_viol_pct(i) = (max_viol / u_range) * 100;

    % --- 5. Estimation Error (% of output magnitude) ---
    % Use c_best(t+1): the filter selected at step t (the active filter
    % whose MPC solution determined simU(:,t)).
    nc = n + r;
    x_hat = zeros(n, steps_sim);
    y_hat = zeros(p, steps_sim);
    for t = 1:steps_sim
        idx = c_best(t+1);
        x_hat(:, t) = RKFs(idx).x_pred(1:n, t);
        y_hat(:, t) = model_nom.C * RKFs(idx).x_pred(1:nc, t);
    end
    err_est = y_hat(:, 1:steps_sim) - simY(:, 1:steps_sim);
    % Normalize by the RMS of the measured output (avoids division by zero at start)
    y_rms = sqrt(mean(simY(:, t_settle:steps_sim).^2, 'all'));
    rmse_est = sqrt(mean(err_est(:, t_settle:end).^2, 'all'));
    estimation_pct(i) = (rmse_est / y_rms) * 100;

    if mod(i, 5) == 0
        fprintf('  Completed %d / %d\n', i, N_mc);
    end
end
sim_time = toc(sim_time);

%% Results Summary
fprintf('\n============================================\n');
fprintf('  PERFORMANCE REPORT  (%d Monte Carlo runs)\n', N_mc);
fprintf('  Reference level = %.1f | Actuator range = [%.1f, %.1f]\n', ...
    ref_level, model_nom.u_min(1), model_nom.u_max(1));
fprintf('Uncertainty level\n  delta_process = %.2f, delta_measure = %.2f\n', ...
    noise_config.delta_process, noise_config.delta_measure);
fprintf('  param_pert_scale = %.2f (element-wise std dev on A, K)\n', param_pert_scale);
fprintf("Controller strategy\n  rkf/lfm = %s\t c-selection = %s\n  offset-free = %s\n", ...
    con_params.mpc, con_params.c_selec, string(MPC_config.offset_free))
fprintf('============================================\n\n');

fprintf('1. TRACKING ERROR (after step at t=%d, %% of reference level = %.1f)\n', t_step, ref_level);
fprintf('   Mean: %5.1f%%  |  Worst-case: %5.1f%%\n', ...
    mean(tracking_pct), max(tracking_pct));
fprintf('   -> "On average, the output deviates %.1f%% from the target"\n\n', ...
    mean(tracking_pct));

fprintf('2. STEADY-STATE ERROR (after settling, t > %d, %% of reference level)\n', t_settle);
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

fprintf('5. ESTIMATION ERROR (after settling, %% of output magnitude)\n');
fprintf('   Mean: %5.1f%%  |  Worst-case: %5.1f%%\n', ...
    mean(estimation_pct), max(estimation_pct));
fprintf('   -> "The filter estimates the output within %.1f%% accuracy"\n\n', ...
    mean(estimation_pct));

fprintf('============================================\n');
fprintf('  Completed in %.2f seconds\n', sim_time);
fprintf('============================================\n');

%% Identify best and worst performing simulations
[~, idx_best]  = min(tracking_pct);
[~, idx_worst] = max(tracking_pct);

fprintf('\n  Best run:  #%2d  (tracking = %5.1f%%, ss = %5.1f%%)\n', ...
    idx_best, tracking_pct(idx_best), ss_error_pct(idx_best));
fprintf('  Worst run: #%2d  (tracking = %5.1f%%, ss = %5.1f%%)\n', ...
    idx_worst, tracking_pct(idx_worst), ss_error_pct(idx_worst));

%% Re-run best and worst simulations for trajectory data
% Deterministic: rng(i) reproduces the exact same noise and perturbation.
fprintf('\nRe-running best and worst for trajectory plots...\n');
traj_best  = struct();
traj_worst = struct();

for ri = [idx_best, idx_worst]
    rng(ri);
    [ms_r, mn_r, ic_r] = model_config(noise_config, MPC_config);
    nr = size(ms_r.A, 1);  mr = size(ms_r.K, 2);
    if param_pert_scale > 0
        for attempt = 1:50
            Ap = ms_r.A .* (1 + param_pert_scale * randn(nr));
            Kp = ms_r.K .* (1 + param_pert_scale * randn(nr, mr));
            if max(abs(eig(Ap))) < 1
                ms_r.A = Ap;  ms_r.K = Kp;
                break;
            end
        end
    end
    [~, sY, tY, sU, ~, ~, ~] = closed_loop_simulation( ...
        ms_r, mn_r, steps_sim, ic_r, reference_ext, set_c, ...
        noise_config.measure_noise, noise_config.process_noise, false, con_params);

    rd = struct( ...
        'simY',  sY(:, 1:steps_sim), ...
        'trueY', tY(:, 1:steps_sim), ...
        'simU',  sU(:, 1:steps_sim), ...
        'track', tracking_pct(ri), ...
        'ss',    ss_error_pct(ri), ...
        'ctrl',  control_pct(ri), ...
        'est',   estimation_pct(ri), ...
        'viol',  constraint_viol_pct(ri));
    if ri == idx_best;  traj_best  = rd;
    else;               traj_worst = rd;
    end
end
fprintf('Done.\n');

%% --- Figure 1: Monte Carlo Performance Distributions ---
col_best  = [0.18 0.62 0.34];
col_worst = [0.84 0.20 0.20];
col_pts   = [0.48 0.48 0.52];
col_box   = [0.82 0.86 0.93];

figure('Name', 'Monte Carlo Performance', ...
    'Position', [50 250 1100 520], 'Color', 'w');
tl1 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, sprintf('Monte Carlo Performance Distributions  (N_{MC} = %d)', N_mc), ...
    'FontSize', 13, 'FontWeight', 'bold');

metric_data   = {tracking_pct, ss_error_pct, control_pct, ...
                 estimation_pct, constraint_viol_pct};
metric_titles = {'Tracking Error', 'Steady-State Error', 'Control Effort', ...
                 'Estimation Error', 'Constraint Violation'};
metric_units  = {'% of reference', '% of reference', '% of actuator range', ...
                 '% of output magnitude', '% of actuator range'};

for k = 1:5
    nexttile; hold on;
    d = metric_data{k};

    boxchart(d, 'BoxFaceColor', col_box, 'MarkerStyle', 'none', 'BoxWidth', 0.45);
    % Individual MC runs with horizontal jitter
    jit = (rand(N_mc, 1) - 0.5) * 0.22;
    scatter(ones(N_mc, 1) + jit, d, 22, col_pts, 'filled', 'MarkerFaceAlpha', 0.55);
    % Highlight best and worst
    scatter(1, d(idx_best),  100, col_best,  'filled', 'pentagram', ...
        'MarkerEdgeColor', col_best * 0.5,  'LineWidth', 0.8);
    scatter(1, d(idx_worst), 100, col_worst, 'filled', 'pentagram', ...
        'MarkerEdgeColor', col_worst * 0.5, 'LineWidth', 0.8);

    hold off; grid on; box on;
    title(metric_titles{k}, 'FontSize', 10);
    ylabel(metric_units{k}, 'FontSize', 9);
    set(gca, 'XTick', [], 'FontSize', 9);
    % Annotate mean and std
    yl = ylim;
    text(1.32, yl(1) + 0.85*(yl(2) - yl(1)), ...
        sprintf('\\mu = %.1f\n\\sigma = %.1f', mean(d), std(d)), ...
        'FontSize', 8, 'Color', [0.35 0.35 0.35]);
end

% Tile 6: legend and configuration summary
nexttile; hold on; axis off;
h1 = scatter(NaN, NaN, 22,  col_pts,   'filled');
h2 = scatter(NaN, NaN, 100, col_best,  'filled', 'pentagram');
h3 = scatter(NaN, NaN, 100, col_worst, 'filled', 'pentagram');
legend([h1 h2 h3], ...
    {'MC run', sprintf('Best  (#%d)', idx_best), sprintf('Worst (#%d)', idx_worst)}, ...
    'Location', 'north', 'FontSize', 10, 'Box', 'off');
text(0.5, 0.22, ...
    {sprintf('Steps = %d', steps_sim), ...
     sprintf('Strategy: %s / %s', con_params.mpc, con_params.c_selec), ...
     sprintf('Offset-free: %s', string(MPC_config.offset_free)), ...
     sprintf('\\delta_p = %.2f   \\delta_m = %.2f', ...
         noise_config.delta_process, noise_config.delta_measure), ...
     sprintf('\\sigma_{param} = %.2f', param_pert_scale)}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Color', [0.45 0.45 0.45]);
hold off;

%% --- Figure 2: Best vs Worst Trajectory Comparison ---
time_ax = 1:steps_sim;
out_colors = lines(max(p, 2));
in_colors  = lines(max(m, 2));

figure('Name', 'Best vs Worst Trajectories', ...
    'Position', [100 60 1100 620], 'Color', 'w');
tl2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Best vs Worst Run — Trajectory Comparison', ...
    'FontSize', 13, 'FontWeight', 'bold');

run_structs = {traj_best, traj_worst};
run_titles  = {sprintf('Best Run (#%d)', idx_best), ...
               sprintf('Worst Run (#%d)', idx_worst)};
run_edge    = {col_best, col_worst};

for col = 1:2
    rd = run_structs{col};

    % --- Row 1: Output tracking ---
    nexttile(col); hold on;
    h_lines = gobjects(p + 1, 1);
    labels  = cell(p + 1, 1);
    for j = 1:p
        % Measured output (faded)
        faded = out_colors(j,:) * 0.3 + [1 1 1] * 0.7;
        plot(time_ax, rd.simY(j,:), ':', 'Color', faded, 'LineWidth', 1);
        % True output (solid)
        h_lines(j) = plot(time_ax, rd.trueY(j,:), '-', ...
            'Color', out_colors(j,:), 'LineWidth', 1.6);
        labels{j} = sprintf('y_%d', j);
    end
    h_lines(end) = plot(time_ax, reference(1, 1:steps_sim), '--k', 'LineWidth', 1.1);
    labels{end} = 'Reference';
    hold off; grid on; box on;
    xline(t_step, ':', 'Color', [0.6 0.6 0.6]);
    xlabel('Time step'); ylabel('Output');
    title(sprintf('%s — Output', run_titles{col}), ...
        'Color', run_edge{col}, 'FontSize', 11);
    legend(h_lines, labels, 'Location', 'best', 'FontSize', 8);
    % Metrics annotation
    yl = ylim;  xl = xlim;
    text(xl(2) - 0.02*(xl(2)-xl(1)), yl(2) - 0.03*(yl(2)-yl(1)), ...
        sprintf('Track: %.1f%%  |  SS: %.1f%%  |  Est: %.1f%%', ...
            rd.track, rd.ss, rd.est), ...
        'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1], 'EdgeColor', [0.75 0.75 0.75], 'Margin', 3);

    % --- Row 2: Control input ---
    nexttile(col + 2); hold on;
    h_ctrl = gobjects(m, 1);
    lbl_ctrl = cell(m, 1);
    for j = 1:m
        h_ctrl(j) = stairs(time_ax, rd.simU(j,:), ...
            'Color', in_colors(j,:), 'LineWidth', 1.4);
        lbl_ctrl{j} = sprintf('u_%d', j);
    end
    % Constraint bounds
    yline(model_nom.u_max(1), '--', 'Color', [0.7 0.15 0.15], ...
        'LineWidth', 1, 'Label', 'u_{max}');
    yline(model_nom.u_min(1), '--', 'Color', [0.7 0.15 0.15], ...
        'LineWidth', 1, 'Label', 'u_{min}');
    hold off; grid on; box on;
    xline(t_step, ':', 'Color', [0.6 0.6 0.6]);
    xlabel('Time step'); ylabel('Control input');
    title(sprintf('%s — Control', run_titles{col}), ...
        'Color', run_edge{col}, 'FontSize', 11);
    legend(h_ctrl, lbl_ctrl, 'Location', 'best', 'FontSize', 8);
    % Metrics annotation
    yl = ylim;  xl = xlim;
    text(xl(2) - 0.02*(xl(2)-xl(1)), yl(2) - 0.03*(yl(2)-yl(1)), ...
        sprintf('Effort: %.1f%%  |  Viol: %.1f%%', rd.ctrl, rd.viol), ...
        'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1], 'EdgeColor', [0.75 0.75 0.75], 'Margin', 3)
end
