function [model_sim, model_nom, init_con] = models(n_model, delta)
arguments
    n_model (1,1) double {mustBeMember(n_model,0:7)}
    delta   (1,1) double {mustBeNonnegative(delta)}
end
switch n_model
    case 0
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [0.5; 1];                     % input -> state
        % model_sim.J = [0.1; 0.05];                % input -> output
    case 1
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [0.5 0; 1 0.1];               % input -> state
        % model_sim.J = [0.1; 0.05];                % input -> output
    case 2
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [0.5 0; 1 0.1];                       % input -> state
        % model_sim.J = [0.1; 0.05];                        % input -> output
    case 3
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [0.5; 1];                             % input -> state
        % model_sim.J = [0.1; 0.05];                        % input -> output
    case 4
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [5; 1];                       % input -> state
        % model_sim.J = [0.1; 0.05];                % input -> output
    case 5
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [5 0.5; 1 0.1];               % input -> state
        % model_sim.J = [0.1; 0.05];                % input -> output
    case 6
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [5 0.5; 1 0.1];                       % input -> state
        % model_sim.J = [0.1; 0.05];                        % input -> output
    case 7
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [5; 1];                       % input -> state
        % model_sim.J = [0.1; 0.05];                        % input -> output
end

% MPC config
model_nom.u_min = -2 * ones(1, size(model_sim.K,2));
model_nom.u_max = 2 * ones(1, size(model_sim.K,2));
model_nom.x_min = [-inf; -inf];
model_nom.x_max = [+inf; +inf];
model_nom.weights.Q = 1*eye(size(model_sim.C,1));
model_nom.weights.Pf = 1*eye(size(model_sim.C,1));
model_nom.weights.R = 0.01*eye(size(model_sim.K,2));

%   init_con:   initial condition
init_con = [1; 0];

%   model_nom:  nominal (perturbed) model used by MPC controller. The
%   perturbation of each entry is the product of the gain delta and a
%   random matrix with compatible sie
model_nom.A = model_sim.A + delta * randn(size(model_sim.A));
model_nom.B = model_sim.B + delta * randn(size(model_sim.B));
model_nom.C = model_sim.C + delta * randn(size(model_sim.C));
model_nom.D = model_sim.D + delta * randn(size(model_sim.D));

model_nom.K = model_sim.K + delta * randn(size(model_sim.K));
% model_con.J = model_sim.J + delta * randn(size(model_sim.J));
end