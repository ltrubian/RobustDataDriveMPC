% PLOTTING Script to plot the results of the closed-loop simulation,
% including control effort, output tracking, and estimation errors.

%% REPORT AND ANALYSIS
% OUTPUT OF THE SIMULATION:
%   simX:       simulated states
%   simY:       simulated output
%   simU:       control input
%   cpuT:       cpu time of the controller
%   filters:    struct with the dynamics of the set of filters
%   c_index:    the sequence of c's selected by the controller
%   x_hat:      state estimation

time_steps = 1:steps_sim;

% Preallocate arrays
x_hat   = zeros(n, steps_sim);
y_hat   = zeros(p, steps_sim);
c_vals  = zeros(1, steps_sim);
sigma_3 = zeros(p, steps_sim);

% Extract the important quantities associated to the selected c-value at each time step
for t = time_steps
    % Identify active filter index (c_index(t+1) is the filter selected at step t)
    idx = c_index(t+1); 
    
    % Extract state prediction and calculate estimated output
    x_hat(:, t) = filters(idx).x_pred(1:n, t);
    y_hat(:, t) = model_nom.C * filters(idx).x_pred(1:n+r, t);
    
    % Extract parameter 'c' and set 0 to 1e-8 for semilogy plotting
    c_vals(t) = max(filters(idx).c, 1e-8);
    
    % Extract 3-Sigma bounds from the least-favorable covariance matrix V
    V_t = model_nom.C * filters(idx).V(1:n+r, 1:n+r, t) * model_nom.C' ...
        + model_nom.D * model_nom.D';
    sigma_3(:, t) = 3 * sqrt(diag(V_t));
end

% Compute state estimation error and output tracking residuals
% NB: the initial state x0 is included in the sequence of true states; 
% the first state in x_true corresponds to the first estimated state in x_hat
x_true = simX(:, 1:steps_sim);
err_x  = x_true - x_hat;
err_track  = reference(:, 1:steps_sim) - simY(:, 1:steps_sim);
err_y = y_hat(:, 1:steps_sim) - simY(:, 1:steps_sim);

%% PLOT FORMATTING SETTINGS
% Centralized settings for consistent aesthetics
lw = 1.5;           % Line width
fs_label = 12;      % Label fontsize
fs_title = 14;      % Title fontsize
latex_opt = {'Interpreter', 'latex'};

% Generate a dynamic color palette large enough for max(n, m, p)
colors = lines(max([n, m, p, 2]));

%% FIGURE 1: Control & Tracking Performance
fig1 = figure('Name', 'System Performance');
tl1 = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1A: Output Tracking
ax1 = nexttile;
hold(ax1, 'on'); 
grid(ax1, 'on');
for i = 1:p
    c = colors(i, :);
    plot(time_steps, reference(i, 1:steps_sim), '--', 'Color', c, 'LineWidth', lw, 'DisplayName', sprintf('Ref $y_%d$', i));
    plot(time_steps, simY(i, 1:steps_sim), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', sprintf('True $y_%d$', i));
    plot(time_steps, y_hat(i, :), ':', 'Color', c, 'LineWidth', lw*1.5, 'DisplayName', sprintf('Est $y_%d$', i));
end
ylabel('Output $y$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{Output Tracking Performance}', latex_opt{:}, 'FontSize', fs_title);
legend('Location', 'best', latex_opt{:});

% 1B: Control Effort
ax2 = nexttile;
hold(ax2, 'on'); 
grid(ax2, 'on');
for i = 1:m
    c = colors(i, :);
    stairs(time_steps, simU(i, 1:steps_sim), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', sprintf('$u_%d$', i));
    if isfield(model_nom, 'u_max') && length(model_nom.u_max) >= i
        yline(model_nom.u_max(i), '--', 'Color', c, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    if isfield(model_nom, 'u_min') && length(model_nom.u_min) >= i
        yline(model_nom.u_min(i), '--', 'Color', c, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
ylabel('Control Input $u$', latex_opt{:}, 'FontSize', fs_label);
xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{Control Effort (with constraints)}', latex_opt{:}, 'FontSize', fs_title);
legend('Location', 'best', latex_opt{:});
linkaxes([ax1, ax2], 'x');
xlim(ax1, [1, steps_sim]);

%% FIGURE 2: Estimation & Adaptive Robustness
fig2 = figure('Name', 'Estimation & Adaptivity');
tl2 = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 2A: State Estimation Error
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on');
for i = 1:n
    plot(time_steps, err_x(i, :), '-', 'Color', colors(i, :), 'LineWidth', lw, 'DisplayName', sprintf('$e_%d$', i));
end
yline(0, 'k--', 'HandleVisibility', 'off');
ylabel('Error $x - \hat{x}$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{State Estimation Error}', latex_opt{:}, 'FontSize', fs_title);
legend('Location', 'best', latex_opt{:});

% 2B: Adaptive Parameter Selection
ax4 = nexttile;
hold(ax4, 'on'); grid(ax4, 'on');
stairs(time_steps, c_vals, 'k-', 'LineWidth', lw);
plot(time_steps, c_vals, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
set(ax4, 'YScale', 'log');
ylabel('Robust Parameter $c$', latex_opt{:}, 'FontSize', fs_label);
xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{Adaptive Selection of Risk-Sensitivity $c$}', latex_opt{:}, 'FontSize', fs_title);
linkaxes([ax3, ax4], 'x');
xlim(ax3, [1, steps_sim]);

%% FIGURE 3: Stochastic Analysis (3-Sigma Bounds)
fig3 = figure('Name', 'Filter Validation: 3-Sigma Bounds');
tl3 = tiledlayout(p, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:p
    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on');
    
    % Plot True Error and 3-Sigma Bounds
    c = colors(i, :);
    plot(time_steps, err_y(i, :), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', 'Estimation Error');
    plot(time_steps, sigma_3(i, :), 'r--', 'LineWidth', 1.2, 'DisplayName', '$\pm 3\sigma$ Bound');
    plot(time_steps, -sigma_3(i, :), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    
    % Statistical calculation for out-of-bounds samples
    out_of_bounds = sum(abs(err_y(i, :)) > sigma_3(i, :));
    perc_out = (out_of_bounds / steps_sim) * 100;
    
    % Dynamic Title & Labels
    title_str = sprintf('\\textbf{Output Estimation Error $\\tilde{y}_%d$ (Out of bounds: %.1f\\%%)}', i, perc_out);
    title(title_str, latex_opt{:}, 'FontSize', fs_title);
    ylabel(sprintf('$\\tilde{y}_%d$', i), latex_opt{:}, 'FontSize', fs_label);
    if i == 1
        legend('Location', 'best', latex_opt{:});
    end
    if i == p
        xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
    end
end
linkaxes(findobj(fig3, 'type', 'axes'), 'x');
xlim(findobj(fig3, 'type', 'axes'), [1, steps_sim]);

%% FIGURE 4: Tracking Error
fig4 = figure('Name', 'Controller Validation: Tracking Error');
tl4 = tiledlayout(p, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:p
    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on');
    
    c = colors(i, :);
    plot(time_steps, err_track(i, :), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', 'Tracking Error: $y_{ref} - y_{sim}$');
    
    yline(0, 'k-', 'HandleVisibility', 'off', 'LineWidth', 0.5);
    
    title_str = sprintf('\\textbf{Tracking Error $y_{ref,%d} - y_%d$}', i, i);
    title(title_str, latex_opt{:}, 'FontSize', fs_title);
    ylabel(sprintf('$e_{track,%d}$', i), latex_opt{:}, 'FontSize', fs_label);
    
    if i == 1
        legend('Location', 'best', latex_opt{:});
    end
    if i == p
        xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
    end
end
linkaxes(findobj(fig4, 'type', 'axes'), 'x');
xlim(findobj(fig4, 'type', 'axes'), [1, steps_sim]);

