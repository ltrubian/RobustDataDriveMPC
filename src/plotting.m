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
sigma_3 = zeros(n, steps_sim);

% Extract the important quantities associated to the selected c-value at each time step
for t = time_steps
    % Identify active filter index
    idx = c_index(t); 
    
    % Extract state prediction and calculate estimated output
    x_hat(:, t) = filters(idx).x_pred(1:n, t);
    y_hat(:, t) = model_nom.C * filters(idx).x_pred(1:n+r, t);
    
    % Extract parameter 'c' and set 0 to 1e-8 for semilogy plotting
    c_vals(t) = max(filters(idx).c, 1e-8);
    
    % Extract 3-Sigma bounds from the least-favorable covariance matrix V
    V_t = filters(idx).V(1:n, 1:n, t);
    sigma_3(:, t) = 3 * sqrt(diag(V_t));
end

% Compute state estimation error and output tracking residuals
% NB: the initial state x0 is included in the sequence of true states; 
% the first state in x_true corresponds to the first estimated state in x_hat
x_true = simX(:, 1:steps_sim);
err_x  = x_true - x_hat;
err_y  = reference(:, 1:steps_sim) - simY(:, 1:steps_sim);

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
linkaxes(findobj(fig3, 'type', 'axes'), 'x');
xlim(findobj(fig3, 'type', 'axes'), [1, steps_sim]);

%% FIGURE 4: MIMO Residual Correlation Analysis
max_lag  = min(20, floor(steps_sim / 4)); 

% Extract full multi-channel matrices [channels x time_steps]
u_data   = simU(:, 1:steps_sim);
eps_data = err_y(:, 1:steps_sim);

% Run MIMO correlation analysis
[R_eps, lag_eps, R_epsu, lag_epsu, conf_limit] = analizza_residui_mimo(eps_data, u_data, max_lag);

% --- FIGURE 4A: Residual Whiteness & Cross-Coupling (p x p Grid) ---
fig4a = figure('Name', 'MIMO Residual Whiteness (p x p)');
tl4a = tiledlayout(p, p, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl4a, '\textbf{MIMO Residual Autocorrelation \& Cross-Coupling $\hat{R}_{\epsilon\epsilon}(\tau)$}', latex_opt{:}, 'FontSize', fs_title);

for i = 1:p
    for j = 1:p
        ax = nexttile;
        hold(ax, 'on'); grid(ax, 'on');
        
        % Extract trajectory across the 3rd dimension
        r_ij = squeeze(R_eps(i, j, :));
        
        stem(ax, lag_eps, r_ij, 'filled', 'Color', colors(i, :), 'LineWidth', lw-0.5, 'MarkerSize', 3);
        yline(ax, conf_limit, 'r--', 'LineWidth', 1);
        yline(ax, -conf_limit, 'r--', 'LineWidth', 1);
        
        ylim(ax, [-1 1]);
        xlim(ax, [-max_lag, max_lag]);
        
        % Labeling matrix coordinates cleanly
        if i == 1
            title(ax, sprintf('Col $\\epsilon_%d$', j), latex_opt{:}, 'FontSize', fs_label);
        end
        if j == 1
            ylabel(ax, sprintf('Row $\\epsilon_%d$', i), latex_opt{:}, 'FontSize', fs_label);
        end
        if i == p
            xlabel(ax, 'Lag $\tau$', latex_opt{:}, 'FontSize', fs_label-2);
        end
    end
end

% --- FIGURE 4B: Residual-Input Independence (p x m Grid) ---
fig4b = figure('Name', 'MIMO Residual-Input Independence (p x m)');
tl4b = tiledlayout(p, m, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl4b, '\textbf{MIMO Residual-Input Cross-Correlation $\hat{R}_{\epsilon u}(\tau)$}', latex_opt{:}, 'FontSize', fs_title);

for i = 1:p
    for j = 1:m
        ax = nexttile;
        hold(ax, 'on'); grid(ax, 'on');
        
        % Extract trajectory across the 3rd dimension
        r_ij = squeeze(R_epsu(i, j, :));
        
        stem(ax, lag_epsu, r_ij, 'filled', 'Color', colors(j, :), 'LineWidth', lw-0.5, 'MarkerSize', 3);
        yline(ax, conf_limit, 'r--', 'LineWidth', 1);
        yline(ax, -conf_limit, 'r--', 'LineWidth', 1);
        
        ylim(ax, [-1 1]);
        xlim(ax, [-max_lag, max_lag]);
        
        % Labeling matrix coordinates cleanly
        if i == 1
            title(ax, sprintf('Input $u_%d$', j), latex_opt{:}, 'FontSize', fs_label);
        end
        if j == 1
            ylabel(ax, sprintf('Res $\\epsilon_%d$', i), latex_opt{:}, 'FontSize', fs_label);
        end
        if i == p
            xlabel(ax, 'Lag $\tau$', latex_opt{:}, 'FontSize', fs_label-2);
        end
    end
end

%% LOCAL FUNCTIONS
function [R_eps, lag_eps, R_epsu, lag_epsu, conf_limit] = analizza_residui_mimo(eps, u, max_lag)
% ANALIZZA_RESIDUI_MIMO Computes normalized auto- and cross-correlation matrices for MIMO systems.
%
% Input:
%   eps     : Residual matrix of size [p, N] (p channels, N time steps)
%   u       : Control input matrix of size [m, N] (m channels, N time steps)
%   max_lag : Maximum lag (\tau) to compute
%
% Output:
%   R_eps      : 3D array [p, p, 2*max_lag+1] of residual auto/cross-correlations
%   lag_eps    : Lag vector for residual correlations
%   R_epsu     : 3D array [p, m, 2*max_lag+1] of residual-input cross-correlations
%   lag_epsu   : Lag vector for input cross-correlations
%   conf_limit : Asymptotic 99% confidence boundary (scalar)

% Validate dimensions
[p, N_eps] = size(eps);
[m, N_u]   = size(u);

if N_eps ~= N_u
    error('Matrices eps and u must have the exact same number of columns (time steps N).');
end
N = N_eps;

% Preallocate 3D correlation arrays
num_lags = 2 * max_lag + 1;
R_eps  = zeros(p, p, num_lags);
R_epsu = zeros(p, m, num_lags);

% 1. Compute p x p Residual Autocorrelation & Cross-Channel Correlation
for i = 1:p
    for j = 1:p
        [r, lag_eps] = xcorr(eps(i, :), eps(j, :), max_lag, 'coeff');
        R_eps(i, j, :) = r;
    end
end

% 2. Compute p x m Residual-Input Cross-Correlation
for i = 1:p
    for j = 1:m
        [r, lag_epsu] = xcorr(eps(i, :), u(j, :), max_lag, 'coeff');
        R_epsu(i, j, :) = r;
    end
end

% 99% asymptotic confidence limit for normalized white noise sequences
conf_limit = 2.58 / sqrt(N);
end