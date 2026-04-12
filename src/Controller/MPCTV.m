function [u_opt, cost_opt] = MPCTV(x0, A, B, C, weights, N, reference, x_min, x_max, u_min, u_max, options)
%MPCOptimizer compute the solution of the MPC problem for linear case
%
%       <usage here>
%
% INPUT:
%   x0          initial condition
%   A           time variant    state->state matrix
%   B           time INvariant  input->state matrix
%   C           time variant    state->output matrix
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
    x0          (:,1)   double
    A           (:,:,:) double
    B           (:,:)   double
    C           (:,:,:) double
    weights     (1,1)   struct
    N           (1,1)   double
    reference   (:,1)   double
    x_min       (:,1)   double
    x_max       (:,1)   double
    u_min       (:,1)   double
    u_max       (:,1)   double
    options
end
arguments(Output)
    u_opt       (:,1) double
    cost_opt    (1,1) double
end

[na, ~, timeA] = size(A);
[~, nc, timeC] = size(C);
[nb, m] = size(B);

assert(na == nb && na == nc, "state dimension do not match")
assert(timeA == N-1 && timeC == N, "notice that you need less A's than C's")

n = na;

% estract blocks of A,C to easily construct the block diagonal matrices
A_sparse = cellfun(@sparse, squeeze(num2cell(A, [1, 2])), 'UniformOutput', false);
C_sparse = cellfun(@sparse, squeeze(num2cell(C, [1, 2])), 'UniformOutput', false);

% generate the block matrices for A and C
A_blk = [sparse(n,N*n)
    blkdiag(A_sparse{:}), sparse((N-1)*n,n)];
C_blk = blkdiag(C_sparse{:});

%% MPC controller setup - sparse formulation
% b) MATRICES to construct EQUALITY CONTRAINT
Geq = [speye(N*n) - A_blk, kron(speye(N), -B)];

% d) BOUNDS ON THE OPTIMIZATION VARIABLE
lb = [repmat(x_min, N, 1); repmat(u_min, N, 1)];
ub = [repmat(x_max, N, 1); repmat(u_max, N, 1)];

fy = C_blk' * blkdiag(kron(speye(N-1), weights.Q), weights.Pf);
Hu = kron(speye(N), weights.R);

H = blkdiag(fy * C_blk, Hu);
f = - [fy * reference; sparse(m*N,1)];

Aeq = sparse(Geq);
beq = x0;

[z_opt, cost_opt, flag, solver_info] = quadprog((H+H')/2, f, [], [], Aeq, beq, lb, ub, [], options);

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