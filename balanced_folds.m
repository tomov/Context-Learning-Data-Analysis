% HACK to get folds such that each fold has 1 irr, 1 mod, and 1 add run
% but randomly shuffled
%
function [foldid, kfolds, c] = balanced_folds(runs, subjs, trials, targets)
    assert(size(targets,1) == numel(runs) * numel(subjs) * numel(trials));
    assert(isequal(numel(runs), 9));
    assert(isequal(numel(subjs), 1));
    % which group of 3 runs we're in (1st 2nd or 3rd)
    iter = [1 * ones(numel(trials) * 3, 1); 2 * ones(numel(trials) * 3, 1); 3 * ones(numel(trials) * 3, 1)];
    % which condition the current run is in
    [~,cond] = max(targets, [], 2);
    foldid = NaN(size(targets,1), 1);
    for c = 1:3 % for each condition
        f = randperm(3); % randomly assign its runs to foldid
        for i = 1:3 % for each iteration of that condition
            assert(sum(iter == i & cond == c) == numel(trials) && sum(isnan(foldid(iter == i & cond == c))) == numel(trials));
            foldid(iter == i & cond == c) = f(i); % assign that iteration of that condition (i.e. the specific run) to a random fold
        end
    end
    assert(sum(foldid == 1) == numel(trials) * 3);
    assert(sum(foldid == 2) == numel(trials) * 3);
    assert(sum(foldid == 3) == numel(trials) * 3);
    assert(sum(isnan(foldid)) == 0);
    kfolds = 3;


    % create a cvpartition object with these folds for compatibility with MATLAB's classifiers
    %
    assert(isequal(numel(subjs), 1));
    c = cvpartition(numel(foldid), 'Kfold', 3);
    c.Impl.indices = foldid;
    c.Impl.TestSize = accumarray(foldid, 1)';
    c.Impl.TrainSize = size(foldid, 1) - c.Impl.TestSize;
end
