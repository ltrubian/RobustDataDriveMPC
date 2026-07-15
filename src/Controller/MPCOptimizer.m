function [u_opt, cost_opt] = MPCOptimizer(x0, A, K, C, weights, N, reference, x_min, x_max, u_min, u_max, options)
%MPCOptimizer compute the solution of the MPC problem for linear case
%
%       <usage here>
%
% INPUT:
%   x0          initial condition
%   A           state matrix
%   K           input-state matrix
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
    x0          (:,1)   double
    A           (:,:,:) double
    K           (:,:)   double
    C           (:,:,:) double
    weights     (1,1) struct
    N           (1,1) double
    reference   (:,1) double
    x_min       (:,1) double
    x_max       (:,1) double
    u_min       (:,1) double
    u_max       (:,1) double
    options
end
arguments(Output)
    u_opt       (:,1) double
    cost_opt    (1,1) double
end
% check and set state dimension
[nb, m] = size(K);
[na, ~, timeA] = size(A);
[~, nc, timeC] = size(C);
assert(na == nb && na == nc, "state dimension do not match")
n = na;

if timeA ~= 1
    % construct time-varying block matrices representing the dynamics of
    % the LFM
    assert(timeA == N-1 && timeC == N, "notice that you need less A's than C's")

    % estract blocks of A,C to easily construct the block diagonal matrices
    A_sparse = cellfun(@sparse, squeeze(num2cell(A, [1, 2])), 'UniformOutput', false);
    C_sparse = cellfun(@sparse, squeeze(num2cell(C, [1, 2])), 'UniformOutput', false);

    % generate the block matrices for A and C
    A_blk = [sparse(n,N*n)
        blkdiag(A_sparse{:}), sparse((N-1)*n,n)];
    C_blk = blkdiag(C_sparse{:});
else
    % construct time-invariant dynamics for the nominal model
    A_blk = kron(diag(ones(N-1,1),-1), A);
    C_blk = kron(speye(N), C);
end

% BOUNDS ON THE OPTIMIZATION VARIABLE
lb = [repmat(x_min, N, 1); repmat(u_min, N, 1)];
ub = [repmat(x_max, N, 1); repmat(u_max, N, 1)];

% MATRICES to construct EQUALITY CONTRAINT (state evolution constraint)
Aeq = sparse([speye(N*n) - A_blk, kron(speye(N), -K)]);
beq = x0;

% COST MATRICES
% construct hessian noticing that the cost are filtered by the matrices C
% since the reference is on the output
fy = C_blk' * blkdiag(kron(speye(N-1), weights.Q), weights.Pf);
Hu = kron(speye(N), weights.R);
% NOTE: H could be used directly, but (H+H')/2 is taken instead to
%       ensure the Hessian matrix to be symmetric even in presence 
%       of numerical errors
f = - [fy * reference; sparse(m*N,1)];
H = blkdiag(fy * C_blk, Hu);
H = (H+H')/2;

idx = 1 + n*N; % starting index for optimal input u 

% selection of the quadratic solver
if ~isempty(options)
    [z_opt, cost_opt, flag, solver_info] = quadprog(H, f, [], [], ...
        Aeq, beq, lb, ub, zeros((n+m)*N,1), options);

    % check the flag to make sure that a solution exists, otherwise, throw error
    if(flag ~= 1)
        error(solver_info.message)
    end
    u_opt  = z_opt(idx:end);
    u_opt = u_opt(1:m);
else
    l = [beq; lb];
    u = [beq; ub];
    A = [Aeq;
        speye(N*(n+m))];

    prob = osqp;
    prob.setup(H, f, A, l, u, 'warm_start', false, 'verbose', false, ...
        'eps_abs', 1e-8, 'eps_rel', 1e-8, 'polish', true);
    res = prob.solve();

    u_opt = res.x(idx:idx+m-1);
    cost_opt = res.info.obj_val;
end
end