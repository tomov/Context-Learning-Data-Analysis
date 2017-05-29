function accuracy = classify_get_accuracy( outputs, targets )
% how well the classifier classified the outputs, compared to the actual
% targets
% each row of targets is a target vector of the form [0 0 ..0 1 0.. 0] 
% each row of outputs is a non-normalized probability distribution, e.g. [1.2 2 3 1 ...]
%

assert(sum(targets(:)) == size(targets, 1));
assert(size(outputs, 1) == size(targets, 1));
assert(size(outputs, 2) == size(targets, 2));

% option 1 -- just see which one is the max
% pretty lame
%
%[~, i] = max(targets, [], 2);
%[~, j] = max(outputs, [], 2);
%accuracy = 100 * mean(i == j);


% option 2 -- average probability of guessing the right answer
%           = 1/n * (P(corr 1) + P(corr 2) + ... + P(corr n)), fundamental bridge =>
%           = 1/n * (E(I1) + E(I2) + ... + E(In)), Ii = indicator r.v. for guessing ith answer correctly
%           = 1/n * E(# right answers)
%           ...this is just an estimator
%
%outputs = outputs ./ sum(outputs, 2); % normalize (just in case) WARNING unsupported by matlab on ncf
assert(abs(mean(sum(outputs, 2)) - 1) < 1e-6);
accuracy = 100 * mean(outputs(logical(targets)));


end

