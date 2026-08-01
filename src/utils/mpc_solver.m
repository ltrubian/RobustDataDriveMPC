function [u_opt, cost_opt] = mpc_solver(x0, A, K, C, weights, N, reference, x_min, x_max, u_min, u_max, options)
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
[p, nc, timeC] = size(C);
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

% MATRICES to construct EQUALITY CONTRAINT (state evolution constraint)
Aeq = sparse([speye(N*n) - A_blk, kron(speye(N), -K)]);
beq = x0;

% Reference computation: compute reference of input u and adapt reference
% of output to the disturabances
BIG_m = [Aeq;C_blk,sparse(p*N,m*N)]; BIG_v = [x0;reference];
ref_xu = BIG_m \ BIG_v;

% COST MATRICES
% construct hessian noticing that the cost are filtered by the matrices C
% since the reference is on the output
fy = C_blk' * blkdiag(kron(speye(N-1), weights.Q), weights.Pf);
Hu = kron(speye(N), weights.R);
% cost on the slack variable must be high so that the solver uses them only
% is strictly necessary in the prediction of the state. The choosen cost is
% 3 order of magnitude higher than the cost on the state
Hs = speye(n*N) * 1e3 * max([weights.Q(:); weights.Pf(:) ] );
% NOTE: H could be used directly, but (H+H')/2 is taken instead to
%       ensure the Hessian matrix to be symmetric even in presence
%       of numerical errors
H = blkdiag(fy * C_blk, Hu, Hs);
f = - H * [ref_xu; ones(N*n, 1)];
H = (H+H')/2;

idx = 1 + n*N; % starting index for optimal input u

% selection of the quadratic solver
if ~isempty(options)
    % BOUNDS ON THE OPTIMIZATION VARIABLE
    lb = [-Inf(N*n, 1); repmat(u_min, N, 1); zeros(N*n, 1)];
    ub = [+Inf(N*n, 1); repmat(u_max, N, 1); +Inf(N*n, 1)];
    Ain = [+speye(N*n), sparse(N*n,N*m), -speye(N*n);
           -speye(N*n), sparse(N*n,N*m), -speye(N*n)];
    bin = [repmat(x_max, N, 1); -repmat(x_min, N, 1)];
    Aeq = [Aeq, sparse(N*n, N*n)];

    [z_opt, cost_opt, flag, solver_info] = quadprog(H, f, Ain, bin, ...
        Aeq, beq, lb, ub, zeros((n+m)*N,1), options);

    % check the flag to make sure that a solution exists, otherwise, throw error
    if(flag ~= 1)
        warning(solver_info.message)
    end
    u_opt  = z_opt(idx:idx+N*m-1);
    % u_opt = u_opt(1:m);
else
    l = [beq; repmat(u_min, N, 1); repmat(x_min, N, 1)];
    u = [beq; repmat(u_max, N, 1); repmat(x_max, N, 1)];
    A = [Aeq, sparse(N*n, N*n);
        sparse(N*m, N*n), speye(N*m), sparse(N*m, N*n);
        speye(N*n), sparse(N*n,N*m), speye(N*n)];

    prob = osqp;
    prob.setup(H, f, A, l, u, 'warm_start', false, 'verbose', false, ...
        'eps_abs', 1e-6, 'eps_rel', 1e-6, 'polish', false);
    res = prob.solve();
    if res.info.status_val ~= 1 && res.info.status_val ~= 2
        warning('OSQP Solver Failed: %s', res.info.status);
    end

    u_opt = res.x(idx:idx+N*m-1);
    cost_opt = res.info.obj_val;
end
end