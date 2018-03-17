function [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, subjs, activations, predict_what, z_score, data, metadata)

% Does the bulk of the work of classify_get_inputs_and_targets
% it's given the betas (activations) pre-computed
%


% condition = context role labels
%
condition_labels = containers.Map(metadata.conditions, {[1 0 0], [0 1 0], [0 0 1]});
                    
n_observations = length(subjs) * length(runs) * length(trials);
n_voxels = size(activations, 2);

% figure out which rows (trials) we're looking at
%
which_rows = ismember(data.participant, metadata.allSubjects(subjs)) & ...
    ismember(data.runId, runs) & ismember(data.newTrialId, trials);
assert(~any(data.drop(which_rows)), 'we should only be using good subjects here; maybe passed some "bad" subjects?');
assert(sum(which_rows) == n_observations, 'maybe wrong order of parameters, e.g. runs and trials');

if strcmp(predict_what, 'responses')
    % if we're looking at subject responses, ignore trials when the subject
    % failed to respond
    which_rows = which_rows & ~data.timeout;
    n_observations = sum(which_rows);
end

%
% Compute input vectors
%

inputs = activations(which_rows, :); % rows = x = observations, cols = voxels / dependent vars
inputs(isnan(inputs)) = 0; % some of them fall outside the imaging range and are NaNs; don't let those screw up our analysis
assert(size(inputs, 1) == n_observations);
assert(size(inputs, 2) == n_voxels);

%
% Compute target vectors
%

targets = []; % rows = y = observations, cols = indep vars (condition as binary vector)

% target depends on what we're trying to predict
%
switch predict_what
    case 'condition'
        targets = values(condition_labels, data.contextRole(which_rows));
        targets = cell2mat(targets);
        
    case 'response'
        targets = zeros(n_observations, 2);
        targets(sub2ind(size(targets), [1:n_observations]', data.chose_sick(which_rows) + 1)) = 1;
        
    case 'runId'
        targets = zeros(n_observations, 9);
        targets(sub2ind(size(targets), [1:n_observations]', data.runId(which_rows))) = 1;
        
    case 'contextId'
        targets = zeros(n_observations, 3);
        targets(sub2ind(size(targets), [1:n_observations]', data.contextId(which_rows) + 1)) = 1;
        
    case 'contextId-training-only'
        assert(all(data.contextId(which_rows) < 2));
        targets = zeros(n_observations, 2);
        targets(sub2ind(size(targets), [1:n_observations]', data.contextId(which_rows) + 1)) = 1;
        
    otherwise
        assert(false, 'should be one of the above');
end            

assert(size(targets, 1) == size(inputs, 1));

%
% Optionally z-score the inputs
%

% get the attributes of the trials we took the betas from
% note that from now on, we'll use these to generate bitmasks referring to
% rows from the the inputs vector, instead of the ones in data
%
runId = data.runId(which_rows);
newTrialId = data.newTrialId(which_rows);
participant = data.participant(which_rows);

switch z_score
    
    case 'z-none'
        % do nothing
        
    case 'z-run'
        inputs = z_run(metadata, runId, newTrialId, participant, subjs, runs, inputs);
        
    case 'z-run-voxel'
        inputs = z_run_voxel(metadata, runId, newTrialId, participant, subjs, runs, trials, inputs);

    case 'pca-subj' 
        inputs = pca_subj(metadata, runId, newTrialId, participant, subjs, runs, inputs);

    case 'z-run-voxel-pca-subj' 
        inputs = z_run_voxel(metadata, runId, newTrialId, participant, subjs, runs, trials, inputs);
        inputs = pca_subj(metadata, runId, newTrialId, participant, subjs, runs, inputs);
        
    otherwise
        assert(false, 'invalid z_score -- should be one of the above');        
end


end % end f'n 



% z-score all voxels within the same run 
% i.e. mean to subtract = mean of all voxels across the run (n_trials_per_run x n_voxels data points)
%
function [inputs] = z_run(metadata, runId, newTrialId, participant, subjs, runs, inputs)
    for subj = subjs
        for run = runs
            which = strcmp(participant, metadata.allSubjects{subj}) & runId == run;
            assert(sum(which) == numel(trials));
            
            inputs_for_run = inputs(which, :);
            z_scored_inputs_for_run = reshape(zscore(inputs_for_run(:)), size(inputs_for_run));
            inputs(which, :) = z_scored_inputs_for_run;
            assert(abs(mean(z_scored_inputs_for_run(:))) < 1e-10);
        end
    end
end

% z-score each voxel within the same run
% i.e. mean to subtract = mean of given voxel across the run (n_trials_per_run data points)
%
function [inputs] = z_run_voxel(metadata, runId, newTrialId, participant, subjs, runs, trials, inputs)
    for subj = subjs
        for run = runs
            which = strcmp(participant, metadata.allSubjects{subj}) & runId == run;
            assert(sum(which) == numel(trials));
            
            inputs_for_run = inputs(which, :);
            inputs(which, :) = zscore(inputs_for_run, 0, 1);
            assert(max(abs(mean(inputs(which, :), 1))) < 1e-10);
        end
    end
end

% PCA for each subject
%
function [inputs] = pca_subj(metadata, runId, newTrialId, participant, subjs, runs, inputs)
    nscores = 10; % how many PCs to include
    for subj = subjs
        which = strcmp(participant, metadata.allSubjects{subj});
        [coeff,score,latent,tsquared,explained,mu] = pca(inputs(which, :));

        new_inputs(which,:) = score(:,1:nscores);
    end
    inputs = new_inputs;
end
