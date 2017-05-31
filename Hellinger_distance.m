function dist = Hellinger_distance( P, Q )
% Compute the Hellinger distance given two PMFs P and Q
% each row is a separate pair of distributions
% https://en.wikipedia.org/wiki/Hellinger_distance
%

assert(isequal(size(P), size(Q)));

dist = 1/sqrt(2) * sqrt(sum((sqrt(P) - sqrt(Q)).^2, 2));

assert(immse(dist, sqrt(1 - exp(-Bhattacharyya_distance(P, Q)))) < 1e-15);


