function [A, B, C, D] = LeastFavorableModel(sys, V, c, N, NF)
%%LeastFavorableModel compute least-favorable model N steps ahead
%
%   [A, B, C, D] = LeastFavorableModel(sys, V, c, N)
%
%   [A, B, C, D] = LeastFavorableModel(sys, V, c, N, NF)
%   if you need to control the length of the forward sweep NF (default= 2N)
%
% the state of returned system A, B, C, D is the extended state composed by
% [x; e], the state x and the estimated error e
%
% WARNING: efficiency is not the goal of this function
%
% INPUT
%   sys:    struct with fields A, B, C, D
%   V:      least-favorable covariance matrix at time 0
%   c:      radius of the ambiguity set
%   N:      how many steps ahead the LFM is computed
%   NF:     how many steps ahead the Riccati iteration on V is computed in
%           order to compute the matrices W, Omega (length forward sweep)
%
% OUTPUT
%   A:      sequence of state -> state matrix
%   B:      sequence of noise -> state matrix
%   C:      sequence of state -> output matrix
%   D:      sequence of noise -> output matrix

arguments
    sys     (1,1) struct
    V       (:,:) double
    c       (1,1) double
    N       (1,1) double {mustBeInteger(N)}
    NF      (1,1) double {mustBeInteger(NF), mustBeGreaterThanOrEqual(NF,N)} = 2*N
end

n = size(sys.A,1);
p = size(sys.C,1);
m = size(sys.B,2);
assert(m == n + p);

% prepare sequence of matrices and update starting point
Vs = zeros(n,n,NF); Vs(:,:,1) = V;
Gs = zeros(n,p,NF);
lambdas = zeros(1,NF);

% forward sweep of risk-sensitive filter
for t=1:(NF - 1)
    [Vs(:,:,t+1), ~, Gs(:,:,t), lambdas(1,t)] = ...
        RiccatiIteration(sys, Vs(:,:,t), c);
end

% backward sweep to evaluate the matrices W, K, H, L
% the inverses of K and Omega are used instead of K and Omega themselves:
% the recursion does not use neither K nor Omega, execept the computation
% of the L matrix where that inversion is done after choleski decomposition
iWs = zeros(n, n, NF);     iWs(:,:,NF) = eye(n)/lambdas(1,NF-1);
iOs = zeros(n, n, NF);
iKs = zeros(n+p, n+p, NF);
Hs = zeros(n+p, n, NF);
Ls = zeros(n+p, n+p, NF);

for t=(NF - 1):-1:1
    % temporary matrices, corrections of sys.A and sys.B
    Bcor = sys.B - Gs(:,:,t) * sys.D;
    Acor = sys.A - Gs(:,:,t) * sys.C;
    % K^-1 and H matrices
    iKs(:,:,t) = (eye(n+p) - Bcor' * iWs(:,:,t+1) * Bcor);
    Hs(:,:,t) = iKs(:,:,t) \ (Bcor' * iWs(:,:,t+1) * Acor);
    % compute L after decomposition of the inverse of K
    Ls(:,:,t) = chol(iKs(:,:,t)) \ eye(size(iKs(:,:,1)));
    % update inverse of Omega to finally compute inverse of K for the next
    % iteration. when t = 1 the cycle is broken because there is no
    % lambda(0) to use (and it is also unnecessary to compute iOs, iWs)
    if t == 1
        break;
    end
    iOs(:,:,t) = Acor' * iWs(:,:,t+1) * Acor + Hs(:,:,t)' * iKs(:,:,t) * Hs(:,:,t);
    iWs(:,:,t) = iOs(:,:,t) + eye(n) / lambdas(1,t-1);
end

% compute the matrices of the extended state
A = zeros(2*n, 2*n, N);
B = zeros(2*n, n+p, N);
C = zeros(  p, 2*n, N);
D = zeros(  p, n+p, N);
for t=1:N
    % temporary matrices, corrections of sys.A and sys.B
    Bcor = sys.B - Gs(:,:,t) * sys.D;
    Acor = sys.A - Gs(:,:,t) * sys.C;
    A(:,:,t) = [
        sys.A,      sys.B * Hs(:,:,t);
        zeros(n),   Acor + Bcor * Hs(:,:,t) ];

    B(:,:,t) = [
        sys.B;
        Bcor ] * Ls(:,:,t);

    C(:,:,t) = [
        sys.C,  sys.D * Hs(:,:,t)];

    D(:,:,t) = sys.D * Ls(:,:,t);
end
end