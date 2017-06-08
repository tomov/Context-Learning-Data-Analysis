function [subjectRDMs, avgSubjectRDM] = compute_rdms(features, distance_measure, data, metadata, which_rows)

% Compute all subject RDMs from a given set of features and distance measure.
%
% INPUT:
% features = [nRows x nFeatures] matrix where nRows = # of rows in data, i.e. these 
%            should correspond exactly to ALL trials in data. We'll be using the logical
%            vectors from data to refer to these rows. These can be, for example, betas
%            for a given mask, or the posteriors from the model, or the stimuli from the task.
% distance_measure = distance argument to pdist, e.g. 'correlation' or
%                    'cosine' or 'euclidean'
% data, metadata = subject behavioral data as output by load_data
% which_rows = which rows to include (make sure that all subjects get the
%              same # of rows)
%

% Compute the single-subject RDMs
%
assert(mod(sum(which_rows), metadata.N) == 0, 'Each subject should have the same number of trials');
trialsPerSubject = sum(which_rows) / metadata.N;
subjectRDMs = nan(trialsPerSubject, trialsPerSubject, metadata.N);

subj_idx = 0;
for who = metadata.subjects
    subj_idx = subj_idx + 1;
    who = who{1};
    
    which_subj_rows = which_rows & strcmp(data.participant, who);
    subjectRDMs(:,:,subj_idx) = squareRDMs(pdist(features(which_subj_rows, :), distance_measure));
end

% Compute the average subject RDM
% display with showRDMs(avgSubjectRDM, 1)
%
avgSubjectRDM = mean(subjectRDMs, 3);

% For displaying with showRDMs(concatSubjectRDMs, 2)
%
%concatSubjectRDMs = cat(3, subjectRDMs, avgSubjectRDM);
