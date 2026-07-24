%% REPORT AND ANALYSIS
% OUTPUT OF THE SIMULATION:
%   simX:       simulated states
%   simY:       simulated output
%   simU:       controlled input
%   cpuT:       cpu time of the controller
%   filters:    struct with the dynamincs of the set of filters
%   c_index:    the sequence of c's selected by the controller
%
% COMPUTED/EXTRACTED VALUES:
%   x_hat:      state estimation

time_steps = 1:steps_sim;

% Preallocate arrays for speed using dynamically extracted dimensions n, p
x_hat   = zeros(n, steps_sim);
y_hat   = zeros(p, steps_sim);
c_vals  = zeros(1, steps_sim);
sigma_3 = zeros(n, steps_sim);

% Extract estimation data, adaptive parameters, and covariances
for t = time_steps
    idx = c_index(t+1); % Identify active filter index

    % Extract state prediction and calculate estimated output
    x_hat(:, t) = filters(idx).x_pred(1:n, t+1);
    y_hat(:, t) = model_con.C(:,1:n+r) * filters(idx).x_pred(1:n+r, t+1);

    % Extract robustness parameter 'c' (set 0 to 1e-8 for semilogy plotting)
    c_vals(t) = max(filters(idx).c, 1e-8);

    % Extract 3-Sigma bounds from the least-favorable covariance matrix V
    V_t = filters(idx).V(1:n, 1:n, t+1);
    sigma_3(:, t) = 3 * sqrt(diag(V_t));
end

% Compute state estimation error (aligning prediction with simulated state)
x_true = simX(:, 2:steps_sim+1);
err_x  = x_true - x_hat;

%% PLOT FORMATTING SETTINGS
% Centralized settings for consistent aesthetics
lw = 1.5;
fs_label = 12;
fs_title = 14;
latex_opt = {'Interpreter', 'latex'};

% Generate a dynamic color palette large enough for max(n, m, p)
colors = lines(max([n, m, p]));

%% 2. FIGURE 1: Control & Tracking Performance
fig1 = figure('Name', 'System Performance');
tl1 = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1A: Output Tracking
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
for i = 1:p
    % Assign consistent color per output
    c = colors(i, :);
    plot(time_steps, reference(i, 1:steps_sim), '--', 'Color', c, 'LineWidth', lw, 'DisplayName', sprintf('Ref $y_%d$', i));
    plot(time_steps, simY(i, 2:steps_sim+1), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', sprintf('True $y_%d$', i));
    plot(time_steps, y_hat(i, :), ':', 'Color', c, 'LineWidth', lw*1.5, 'DisplayName', sprintf('Est $y_%d$', i));
end
ylabel('Output $y$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{Output Tracking Performance}', latex_opt{:}, 'FontSize', fs_title);
legend('Location', 'best', latex_opt{:});

% 1B: Control Effort
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
for i = 1:m
    c = colors(i, :);
    stairs(time_steps, simU(i, 1:steps_sim), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', sprintf('$u_%d$', i));

    % Dynamically plot input constraints if they exist for this dimension
    if isfield(model_con, 'u_max') && length(model_con.u_max) >= i
        yline(model_con.u_max(i), '--', 'Color', c, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    if isfield(model_con, 'u_min') && length(model_con.u_min) >= i
        yline(model_con.u_min(i), '--', 'Color', c, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
ylabel('Control Input $u$', latex_opt{:}, 'FontSize', fs_label);
xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{Control Effort (with constraints)}', latex_opt{:}, 'FontSize', fs_title);
legend('Location', 'best', latex_opt{:});

linkaxes([ax1, ax2], 'x');
xlim(ax1, [1, steps_sim]);

%% 3. FIGURE 2: Estimation & Adaptive Robustness
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

%% 4. FIGURE 3: Stochastic Analysis (3-Sigma Bounds)
% Scale figure height dynamically based on the number of states 'n'
fig3 = figure('Name', 'Filter Validation: 3-Sigma Bounds');
tl3 = tiledlayout(n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:n
    ax = nexttile;
    hold(ax, 'on'); grid(ax, 'on');

    % Plot True Error and 3-Sigma Bounds
    c = colors(i, :);
    plot(time_steps, err_x(i, :), '-', 'Color', c, 'LineWidth', lw, 'DisplayName', 'Estimation Error');
    plot(time_steps, sigma_3(i, :), 'r--', 'LineWidth', 1.2, 'DisplayName', '$\pm 3\sigma$ Bound');
    plot(time_steps, -sigma_3(i, :), 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    % Statistical calculation for out-of-bounds samples
    out_of_bounds = sum(abs(err_x(i, :)) > sigma_3(i, :));
    perc_out = (out_of_bounds / steps_sim) * 100;

    % Dynamic Title & Labels
    title_str = sprintf('\\textbf{State Error $e_%d$ (Out of bounds: %.1f\\%%)}', i, perc_out);
    title(title_str, latex_opt{:}, 'FontSize', fs_title);
    ylabel(sprintf('$e_%d$', i), latex_opt{:}, 'FontSize', fs_label);

    if i == 1
        legend('Location', 'best', latex_opt{:});
    end
    if i == n
        xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
    end
end

% Link all X-axes in Figure 3
linkaxes(findobj(fig3, 'type', 'axes'), 'x');
xlim(findobj(fig3, 'type', 'axes'), [1, steps_sim]);

%%%% RESIDUAL ANALYSIS DIAGNOSTICS (Commented by default; uncomment to run)
%
% % Generate dummy test data for analysis
% u_data = simU(1:end-1);
% eps_data = reference(1,1:steps_sim-1) - simY(1,2:end-1);
% 
% % Function call
% max_lag = 20;
% [r_eps, lag_eps, r_epsu, lag_epsu, conf_limit] = analyze_residuals_normalized(eps_data, u_data, max_lag);
% 
% % Plot Normalized Autocorrelation
% figure;
% subplot(2,1,1);
% stem(lag_eps, r_eps, 'filled');
% hold on;
% yline(conf_limit, '--r', '99% Conf Limit');
% yline(-conf_limit, '--r');
% title('Normalized Autocorrelation of Residuals');
% xlabel('Lag (\tau)');
% ylim([-1 1]); % Set Y limits since the signal is normalized
% grid on;
% 
% % Plot Normalized Cross-correlation
% subplot(2,1,2);
% stem(lag_epsu, r_epsu, 'filled');
% hold on;
% yline(conf_limit, '--r', '99% Conf Limit');
% yline(-conf_limit, '--r');
% title('Normalized Cross-correlation of Residuals and Inputs');
% xlabel('Lag (\tau)');
% ylim([-1 1]);
% grid on;
%
% function [r_eps, lag_eps, r_epsu, lag_epsu, conf_limit] = analyze_residuals_normalized(eps, u, max_lag)
% % ANALYZE_RESIDUALS_NORMALIZED Calculates normalized auto- and cross-correlation.
% %
% % Input:
% %   eps     : Row vector of residuals (filter innovation)
% %   u       : Row vector of control inputs
% %   max_lag : Maximum lag to calculate
% %
% % Output:
% %   r_eps      : Autocorrelation of residuals (normalized, max 1)
% %   lag_eps    : Lag vector for autocorrelation
% %   r_epsu     : Cross-correlation of residuals and inputs (normalized)
% %   lag_epsu   : Lag vector for cross-correlation
% %   conf_limit : Asymptotic 99% confidence limit (scalar)
% 
% % Ensure inputs are row vectors
% eps = eps(:).';
% u = u(:).';
% 
% N = length(eps);
% if length(u) ~= N
%     error('Vectors eps and u must have the same length.');
% end
% 
% % 1. Calculate normalized Autocorrelation
% % 'coeff' normalizes the sequence so that autocorrelation at lag 0 is 1
% [r_eps, lag_eps] = xcorr(eps, max_lag, 'coeff');
% 
% % 2. Calculate normalized Cross-correlation
% % 'coeff' divides cross-correlation by sqrt(var(eps) * var(u))
% [r_epsu, lag_epsu] = xcorr(eps, u, max_lag, 'coeff');
% 
% % Asymptotic 99% confidence limit for normalized signals
% conf_limit = 2.58 / sqrt(N);
% end
