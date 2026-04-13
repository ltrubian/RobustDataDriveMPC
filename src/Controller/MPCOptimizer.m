function [u_opt, cost_opt] = MPCOptimizer(x0, A, B, C, weights, N, reference, x_min, x_max, u_min, u_max)
%MPCOptimizer compute the solution of the MPC problem for linear case
%
%       <usage here>
%
% INPUT:
%   x0          initial condition
%   A           state matrix
%   B           input-state matrix
%   C           state-output matrix
%   weights     weights of the costs (Q,Pf for output, R for input)
%   N           prediction horizon
%   reference   output reference
%   x_min       min value for each state
%   x_max       max value for each state
%   u_min       min value for each input
%   u_max       max value for each input
%
% OUTPUT:
%   u_opt:       input of the system
%   cost_opt:    value of the cost function at its optimum

arguments(Input)
    x0          (:,1) double
    A           (:,:) double
    B           (:,:) double
    C           (:,:) double
    weights     (1,1) struct
    N           (1,1) double
    reference   (:,1) double
    x_min       (:,1) double
    x_max       (:,1) double
    u_min       (:,1) double
    u_max       (:,1) double
end
arguments(Output)
    u_opt       (:,1) double
    cost_opt    (1,1) double
end

n = size(A,1);  % number of state
m = size(B,2);  % number of inputs

%% MPC controller setup - sparse formulation

% a) OBJECTIVE HESSIAN
% H_bar = blkdiag(kron(eye(N-1), C' * weights.Q * C), C' * weights.Pf * C, kron(eye(N), weights.R));

% b) MATRICES to construct EQUALITY CONTRAINT
Geq = [kron(eye(N), eye(n)) - kron(diag(ones(N-1,1),-1), A), kron(eye(N), -B)];
Eeq = [A; zeros((N-1)*n, n)];

% c) MATRICES and VECTOR to construct GENERAL INEQUALITY CONSTRAINT
% In this case this constraint are enforced only on starting x0
Gin = [zeros(2*n, (n+m)*N)];
win = [x_max; -x_min];
Ein = [-eye(n); eye(n)];

% d) BOUNDS ON THE OPTIMIZATION VARIABLE
lb = [repmat(x_min, N, 1); repmat(u_min, N, 1)];
ub = [repmat(x_max, N, 1); repmat(u_max, N, 1)];

H  = sparse(blkdiag(kron(eye(N-1), C' * weights.Q * C), C' * weights.Pf * C, kron(eye(N), weights.R)));
f  = -blkdiag(kron(eye(N-1), C' * weights.Q), C' * weights.Pf, kron(eye(N), weights.R))*[reference; zeros(m*N,1)];

Aeq = sparse(Geq);
beq = Eeq*x0;

Ain = sparse(Gin);
bin = Ein*x0 + win;

% compute optimal input sequence
options = optimset('quadprog');
options = optimset(options, 'Algorithm', 'interior-point-convex', 'Display', 'off');

[z_opt, cost_opt, flag, solver_info] = quadprog((H+H')/2, f, Ain, bin, Aeq, beq, lb, ub, [], options);

%NOTE: H could be used directly, but (H+H')/2 is taken instead to
%      ensure the Hessian matrix to be symmetric even in presence of numerical errors

% check the flag to make sure that a solution exists, otherwise, throw error
if(flag ~= 1)
    error(solver_info.message)
end

idx = 1 + n*N;
u_opt  = z_opt(idx:end);


% apply the first control input sample
u_opt = u_opt(1:m);


end