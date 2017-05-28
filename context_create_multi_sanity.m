function context_create_multi_sanity(glmodels, subjs, runs)

% Sanity test context_create_multi.m before shipping to NCF
%
% INPUT:
% glms = glmodels to test, e.g. 1:20
% subjs (optional) = subject indices to test, e.g. getGoodSubjects()
% runs (optional) = runs to test, e.g. 1:9
%

% set default parameters
%
if nargin < 3
    runs = 1:9;
end

if nargin < 2
    subjs = getGoodSubjects();
end

for glmodel = glmodels
    for subj = subjs
        for run = runs
            multi = context_create_multi(glmodel, subj, run);
        end
    end
end

