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
%   x_hat:      state estimation made by the (combination of) controller and
%               filters and used by the first to produce the input and the
%               latter estimate the next state

% the struct filters holds all the prediction made by all the filters and
% we need to extract the ones that were actually used by the controller

time_steps = 1:steps_sim;

% Preallocate arrays for speed
x_hat   = zeros(n, steps_sim);
y_hat   = zeros(p, steps_sim);
c_vals  = zeros(1, steps_sim);
sigma_3 = zeros(n, steps_sim);

% Extract estimation data, adaptive parameters, and covariances
for t = time_steps
    idx = c_index(t+1); % Identify active filter index
    
    % Extract state prediction and calculate estimated output
    x_hat(:, t) = filters(idx).x_pred(1:n, t+1);
    y_hat(:, t) = model_con.C * x_hat(:, t);
    
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

%% 2. FIGURE 1: Control & Tracking Performance
fig1 = figure('Name', 'System Performance');
tl1 = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1A: Output Tracking
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
plot(time_steps, reference(1, 1:steps_sim), 'k--', 'LineWidth', lw);
plot(time_steps, simY(1, 2:steps_sim+1), 'b-', 'LineWidth', lw);
plot(time_steps, y_hat(1, :), 'r-.', 'LineWidth', lw);
ylabel('Output $y$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{Output Tracking Performance}', latex_opt{:}, 'FontSize', fs_title);
legend('Reference', 'True Output', 'Estimated Output', 'Location', 'best', latex_opt{:});

% 1B: Control Effort
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
stairs(time_steps, simU(1, 1:steps_sim), 'b-', 'LineWidth', lw);
stairs(time_steps, simU(2, 1:steps_sim), 'r-', 'LineWidth', lw);
yline(model_con.u_max(1), 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(model_con.u_min(1), 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
ylabel('Control Input $u$', latex_opt{:}, 'FontSize', fs_label);
xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{Control Effort (with constraints)}', latex_opt{:}, 'FontSize', fs_title);
legend('$u_1$', '$u_2$', 'Location', 'best', latex_opt{:});

linkaxes([ax1, ax2], 'x');
xlim(ax1, [1, steps_sim]);

%% 3. FIGURE 2: Estimation & Adaptive Robustness
fig2 = figure('Name', 'Estimation & Adaptivity');
tl2 = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 2A: State Estimation Error
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on');
plot(time_steps, err_x(1, :), 'b-', 'LineWidth', lw);
plot(time_steps, err_x(2, :), 'r-', 'LineWidth', lw);
yline(0, 'k--', 'HandleVisibility', 'off');
ylabel('Error $x - \hat{x}$', latex_opt{:}, 'FontSize', fs_label);
title('\textbf{State Estimation Error}', latex_opt{:}, 'FontSize', fs_title);
legend('$e_1$', '$e_2$', 'Location', 'best', latex_opt{:});

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
    plot(time_steps, err_x(i, :), 'b-', 'LineWidth', lw);
    plot(time_steps, sigma_3(i, :), 'r--', 'LineWidth', 1.2);
    plot(time_steps, -sigma_3(i, :), 'r--', 'LineWidth', 1.2);
    
    % Statistical calculation for out-of-bounds samples
    out_of_bounds = sum(abs(err_x(i, :)) > sigma_3(i, :));
    perc_out = (out_of_bounds / steps_sim) * 100;
    
    % Dynamic Title & Labels
    title_str = sprintf('\\textbf{State Error $e_%d$ (Out of bounds: %.1f\\%%)}', i, perc_out);
    title(title_str, latex_opt{:}, 'FontSize', fs_title);
    ylabel(sprintf('$e_%d$', i), latex_opt{:}, 'FontSize', fs_label);
    
    if i == 1
        legend('Estimation Error', '$\pm 3\sigma$ Bound', 'Location', 'best', latex_opt{:});
    end
    if i == n
        xlabel('Time Step $k$', latex_opt{:}, 'FontSize', fs_label);
    end
end

% Link all X-axes in Figure 3
linkaxes(findobj(fig3, 'type', 'axes'), 'x');
xlim(findobj(fig3, 'type', 'axes'), [1, steps_sim]);