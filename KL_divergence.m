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
logs = log2(P) - log2(Q); 
logs(isnan(logs)) = 0; % lim_{x->0} x log(x) = 0
surprise = sum(P .* logs, 2);
surprise(isnan(surprise)) = 0; % weird things happen when P --> 0, e.g. we get -Infs

end

