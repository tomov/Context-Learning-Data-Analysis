function surprise = KL_divergence(P, Q)
% Compute the KL divergence from Q to P i.e. D_KL(P || Q)
% = sum over i, 
% assumes each column corresponds to a different value of the RV
% if given matrices, treats each row as a separate pair of P, Q
%
% INPUT:
% P = posterior distributions, one per row
% Q = prior distributions, one per row
%
% OUTPUT:
% surprise = D_KLs, one per row
%

assert(isequal(size(P), size(Q)));

h = P .* (log2(P) - log2(Q)); 
h(isnan(h)) = 0; % lim_{x->0} x log(x) = 0
assert(sum(sum(isinf(h))) == 0); % shouldn't be any inf's left

surprise = sum(h, 2);

assert(sum(isinf(surprise)) == 0); % sanity
assert(sum(isnan(surprise)) == 0); % sanity

end

