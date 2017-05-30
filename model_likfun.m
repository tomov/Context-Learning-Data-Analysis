function total_loglik = model_likfun(data, metadata, params, which_structures, which_rows, testTrialsOnly)
% Log likelihood function of the subject data according to the model predictions
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% params = vector of hyperparameters:
%    params(:, 1) = prior_variance
%    params(:, 2) = inv_softmax_temp
%          for fixed effects, only have 1 row and these parameters will be
%          used for all subjects. For random effects, have a separate row
%          with the parameters for each subject
% which_structures = which causal structures to use, as a logical vector,
%                    e.g. [1 1 1 0] = M1, M2, M3
% which_rows = binary mask saying which rows from the data to use. For
%              fixed effects, this will include all subjects. For random
%              effects, it will include one subject at a time.


% First simulate the subjects with the causal structure model
%
simulated = simulate_subjects(data, metadata, params, which_structures, which_rows);

% Exclude timeouts & training trials (optional)
%
which = which_rows & ~data.timeout;
if testTrialsOnly
    which = which & ~data.isTrain;
end

X = data.chose_sick(which); % actual subject choice on each trial
P = simulated.pred(which); % probability of subject choice on each trial

% lik = binopdf(X, 1, P); <-- unsupported by old MATLAB on NCF
% assert(numel(lik) == numel(X));
% total_loglik = sum(log(lik));

total_loglik = 0;
assert(numel(X) == numel(P));
for i = 1:numel(X)
    lik = binopdf(X(i), 1, P(i));
    total_loglik = total_loglik + log(lik);
end