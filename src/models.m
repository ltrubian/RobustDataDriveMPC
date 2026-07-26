function [model_sim, model_nom, init_con] = models(n_model, delta_process, delta_measure, offset_free)
arguments
    n_model       (1,1) double {mustBeMember(n_model,0:7)}
    delta_process (1,1) double {mustBeNonnegative(delta_process)}
    delta_measure (1,1) double {mustBeNonnegative(delta_measure)}
    offset_free   (1,1) logical
end
switch n_model
    case 0
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [0.5; 1];                     % input -> state
    case 1
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [0.5 0; 1 0.1];               % input -> state
    case 2
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [0.5 3; 1 0.1];                       % input -> state
    case 3
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [0.5; 1];                             % input -> state
    case 4
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [5; 1];                       % input -> state
    case 5
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                 % state -> state
        model_sim.B = [0.5 0.2 0.1; 0.3 0.2 0.01];  % noise -> state
        model_sim.C = [1 0];                        % state -> output
        model_sim.D = [0.1, 0.05, 0.01];            % noise -> output
        model_sim.K = [5 0.5; 1 0.1];               % input -> state
    case 6
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [5 0.5; 1 0.1];                       % input -> state
    case 7
        %   model_sim:  true model to simulate
        model_sim.A = [1.1 1; 0 1];                         % state -> state
        model_sim.B = [0.5 0.2 0.1 0; 0 0.3 0.2 0.01];      % noise -> state
        model_sim.C = [1 0; 0.3 1];                         % state -> output
        model_sim.D = [0.1, 0.05, 0.01 0.2; 0 0.2 0.1 0.05];% noise -> output
        model_sim.K = [5; 1];                               % input -> state
end
n = size(model_sim.A, 1);
m = size(model_sim.K, 2);
p = size(model_sim.C, 1);

% MPC config

% bounds
model_nom.u_min = -2 * ones(1, size(model_sim.K,2));
model_nom.u_max = 2 * ones(1, size(model_sim.K,2));
model_nom.x_min = [-inf; -inf];
model_nom.x_max = [+inf; +inf];

% weights
model_nom.weights.Q = 1*eye(size(model_sim.C,1));
model_nom.weights.Pf = 1*eye(size(model_sim.C,1));
model_nom.weights.R = 0.1*eye(size(model_sim.K,2));

%   init_con:   initial condition
init_con = [1; 0];

% model_nom: matrices A,C,K are the same as the real model, matrices B,D
% are respectively the process and measure gain matrices
model_nom.A = model_sim.A;
model_nom.B = [delta*eye(n), zeros(n, p)];
model_nom.C = model_sim.C;
model_nom.D = [zeros(p, n), delta*eye(p)];
model_nom.K = model_sim.K;

% disturbances matrices for the offset-free tracking. Due to detectability
% and the proprerty of offset-free tracking, the number of disturbances is
% equal to the number of tracked output (in our case, we try to track all
% the outputs
tmp_disturbances = (ones(n+p, p) + eye(n+p,p))/2;
model_nom.Bd = tmp_disturbances(1:n,:);
model_nom.Cd = tmp_disturbances(n+1:end,:) + eye(p);

% assert detectability of disturbances
assert(rank([model_nom.A - eye(n), model_nom.Bd; ...
    model_nom.C, model_nom.Cd]) == n+p, "Non-observable disturbances");
% assert dimensions
assert(size(model_sim.K, 1) == n, "invalid #rows input-state matrix");
assert(size(model_sim.C, 2) == n, "invalid #cols state-output matrix");
assert(size(model_sim.B, 1) == n, "invalid #rows noise-state matrix");
assert(size(model_sim.D, 1) == p, "invalid #rows noise-output matrix");
assert(size(model_sim.B, 2) == n+p, "invalid #cols noise-state matrix");
assert(size(model_sim.D, 2) == n+p, "invalid #cols noise-output matrix");
% assert controllability
Co = ctrb(model_sim.A, model_sim.K);
assert(rank(Co) == n, "Non-controllable real world model");
Co = ctrb(model_nom.A, model_nom.K);
assert(rank(Co) == n, "Non-controllable nominal model");
% assert observability
Ob = obsv(model_sim.A, model_sim.C);
assert(rank(Ob) == n, "Non-observable real world model");
Ob = obsv(model_nom.A, model_nom.C);
assert(rank(Ob) == n, "Non-observable nominal model");
% assert noise covariance invertibility
assert(rank([model_sim.B; model_sim.D]) == n+p, "Non-Invertible noise covariance")
assert(rank([model_nom.B; model_nom.D]) == n+p, "Non-Invertible noise covariance")

if offset_free
    model_nom.A = [model_nom.A model_nom.Bd; zeros(p, n) eye(p)];
    model_nom.K = [model_nom.K; zeros(p, m)];
    model_nom.C = [model_nom.C model_nom.Cd];
    model_nom.B = [model_nom.B, zeros(n, p); zeros(p, n+p), eye(p)];
    model_nom.D = [model_nom.D, zeros(p, p)];
    model_nom.x_min = [model_nom.x_min; -Inf(p,1)];
    model_nom.x_max = [model_nom.x_max; +Inf(p,1)];
    assert(rank([model_nom.B; model_nom.D]) == n+2*p, "Non-Invertible noise covariance")
end
end