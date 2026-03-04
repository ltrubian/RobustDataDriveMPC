function A = random_cov(n)
% https://math.stackexchange.com/questions/357980/how-to-generate-random-symmetric-positive-definite-matrices-using-matlab#358092

Q = randn(n,n);

eigen_mean = 1;
% can be made anything, even zero
% used to shift the mode of the distribution

A = Q' * diag(abs(eigen_mean+randn(n,1))) * Q;

return