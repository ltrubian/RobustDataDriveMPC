function [A, B, C, D] = LeastFavorableModel(sys, V, P, N, c)
%%LeastFavorableModel
%
%   [A, B, C, D] = LeastFavorableModel(sys, V, N, c)
%   compute least-favorable model N steps ahead
%
% NOTE: the state of returned system A,B,C,D is the extended state composed
%       by [x; e], the state x and the estimated error e
%
% INPUT
%   sys:    struct with fields A, B, C, D
%   V:      least-favorable covariance matrix at time 0
%   N:      how many steps ahead the LFM is computed
%   c:      radius of the ambiguity set
%
% OUTPUT
%   A:      sequence of state -> state matrix
%   B:      sequence of noise -> state matrix
%   C:      sequence of state -> output matrix
%   D:      sequence of noise -> output matrix

arguments
    sys     (1,1) struct
    V       (:,:) double
    P       (:,:) double
    N       (1,1) double
    c       (1,1) double
end

n = size(sys.A,1);
p = size(sys.C,2);

Vs = zeros(n,n,2*N); Vs(:,:,1) = V;
Ps = zeros(n,n,2*N); Ps(:,:,1) = P;
Gs = zeros(n,p,2*N);
lambdas = zeros(1,2*N);

% forward sweep of risk-sensitive filter
for t=1:(2*N - 1)
    [Vs(:,:,t+1), Ps(:,:,t+1), Gs(:,:,t), lambdas(1,t)] = ...
        RiccatiIteration(sys, V(:,:,t), c);
end
[~,~, Gs(:,:,2*N), lamddas(1,2*N)] = RiccatiIteration(sys, V(:,:,2*N), c);

% TODO: backward sweep


end