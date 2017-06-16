function surprise = KL_divergence_gauss(Mu0, Sigma0, Mu1, Sigma1)
% Compute the KL divergence between two multivariate gaussians
% with means m0 and Mu1 and covariances sigma0 and sigma1
% It's D_KL(gaussian 0 || gaussian 1).
% Does this for N pairs of distributions
%
% INPUTS:
% Mu0, Mu1 = N x D arrays with the means, where N = # of gaussians, D = # variables
% sigma0, sigma1 = D x D x N arrays with the covariance matrices

N = size(Mu0, 1);
D = size(Mu1, 2);

assert(N == size(Mu1, 1));
assert(N == size(Sigma0, 3));
assert(N == size(Sigma1, 3));

assert(D == size(Mu1, 2));
assert(D == size(Sigma0, 1));
assert(D == size(Sigma0, 2));
assert(D == size(Sigma1, 1));
assert(D == size(Sigma1, 2));

surprise = nan(N, 1);
for i = 1:size(Mu0, 1)
    mu0 = Mu0(i,:)';
    mu1 = Mu1(i,:)';
    sigma0 = Sigma0(:,:,i);
    sigma1 = Sigma1(:,:,i);

    surprise(i) = trace(sigma1 \ sigma0) + (mu1 - mu0)' * (sigma1 \ (mu1 - mu0)) - D + log(det(sigma1) / det(sigma0));
    surprise(i) = 0.5 * surprise(i);
    surprise(i) = surprise(i) / log(2); % convert from nats to bits
end

