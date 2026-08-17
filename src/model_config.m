function [model_sim, model_nom, init_con] = model_config(noise_config, MPC_config)
% MODEL_CONFIG Configures simulation and nominal state-space models for MPC.
%
% This function loads the quadruple tank base model and prepares:
%   1. model_sim: The "real" plant model used for simulation.
%   2. model_nom: The internal MPC model with intentional plant-model mismatches,
%                 assigned constraints/weights, and optional disturbance state 
%                 augmentation for offset-free tracking.
%   3. init_con:  The initial state conditions.
% 
% The function also asserts system dimension validity, detectability of 
% disturbances, controllability, observability, and invertibility of noise covariance.

arguments
    noise_config (1,1) struct
    MPC_config   (1,1) struct
end

% unpack noise configuration
delta_process = noise_config.delta_process;
delta_measure = noise_config.delta_measure;
mustBeNonnegative(delta_process);
mustBeNonnegative(delta_measure);

% unpack MPC configuration
offset_free = MPC_config.offset_free;

% generate .mat file with system matrices if not already present
fileName = 'quadruple_tank.mat';
if ~isfile(fileName)
    generate_quadruple_tank_model();
end
load(fileName, 'A', 'B', 'C', 'D', 'K');
model_sim.A = A;
model_sim.B = B;
model_sim.C = C;
model_sim.D = D;
model_sim.K = K;
clear A B C D K

n = size(model_sim.A, 1);
m = size(model_sim.K, 2);
p = size(model_sim.C, 1);

% set initial conditions
init_con = [10; 10; 1; 1];

% matrices A, C, and K include intentional mismatches to demonstrate the 
% effectiveness of offset-free tracking. Matrices B and D define the 
% process and measurement noise gains, respectively.
model_nom.A = model_sim.A * 0.95;
model_nom.B = [delta_process*eye(n), zeros(n, p)];
model_nom.C = model_sim.C * 1.05;
model_nom.D = [zeros(p, n), delta_measure*eye(p)];
model_nom.K = model_sim.K * 1.1;

% build MPC bounds and weights from scalar config in MPC_config
model_nom.u_min = MPC_config.u_min * ones(m, 1);
model_nom.u_max = MPC_config.u_max * ones(m, 1);
model_nom.x_min = MPC_config.x_min * ones(n, 1);
model_nom.x_max = MPC_config.x_max * ones(n, 1);
model_nom.weights.Q  = MPC_config.Q * eye(p);
model_nom.weights.Pf = MPC_config.Pf * eye(p);
model_nom.weights.R  = MPC_config.R * eye(m);

% disturbances matrices for the offset-free tracking. Due to detectability
% and the property of offset-free tracking, the number of disturbances is
% equal to the number of tracked output (in our case, we try to track all
% the outputs)
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

% expand the matrices to allow for offset free tracking
if offset_free
    model_nom.A = [model_nom.A model_nom.Bd; zeros(p, n) eye(p)];
    model_nom.K = [model_nom.K; zeros(p, m)];
    model_nom.C = [model_nom.C model_nom.Cd];
    model_nom.B = [model_nom.B, zeros(n, p); zeros(p, n+p), eye(p)];
    model_nom.D = [model_nom.D, eye(p, p)];
    model_nom.x_min = [model_nom.x_min; -Inf(p,1)];
    model_nom.x_max = [model_nom.x_max; +Inf(p,1)];
    assert(rank([model_nom.B; model_nom.D]) == n+2*p, "Non-Invertible noise covariance")
end
end