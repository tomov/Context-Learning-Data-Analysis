function dist = Bhattacharyya_distance( P, Q )
% Compute the Bhattacharyya distance given two PMFs P and Q
% each row is a separate pair of distributions
% https://en.wikipedia.org/wiki/Bhattacharyya_distance
%

assert(isequal(size(P), size(Q)));

BC = sum(sqrt(P.*Q), 2);
dist = - log(BC);

%assert(immse(Hellinger_distance(P, Q), sqrt(1 - BC)) < 1e-15);

