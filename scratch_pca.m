% PCA stuff based on https://github.com/tomov/neurolab/tree/master/exercise

%% load behavioral stuff

% Load data
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% Load parameters
%
load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts_5_prior.mat'), 'results', 'results_options');
params = results(1).x;
options = results_options(1);

% OVERRIDE -- use pilot params as before
%
params = [0.1249 2.0064];
options.isFmriData = false;
options.fixedEffects = true;

disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);
% safeguards
%assert(options.isFmriData == true);
%assert(~options.fixedEffects);
assert(isequal(options.which_structures, [1 1 1 0]));
which_structures = logical([1 1 1 0]);

% Run the model with the parameters
%
simulated = simulate_subjects(data, metadata, params, which_structures);

%% load betas -- SLOW
%

%whole_brain_betas = get_betas('masks/mask.nii', 'trial_onset', data, metadata, false);

%% load mask & corresponding betas
%
%bspmview('masks/prior_left_IFG.nii', 'masks/mean.nii');

[m, V] = load_mask('masks/prior_left_IFG.nii'); % left IFG for now TODO right one too
betas = get_activations_submask(m, whole_brain_betas);
assert(size(betas, 1) == size(data.which_rows, 1));

%% restrict to relevant trials
%
which_trials = data.which_rows; % & data.isTrain; % Look at training trials only

figure;
title('Data in voxel space');
imagesc(betas(which_trials, :));
xlabel('voxel');
ylable('trial');

%% run PCA 
%
[coeff,score,latent,tsquared,explained,mu] = pca(betas(which_trials, :));
score(which_trials, :) = score; % include dummy trials for easier indexing

figure;
title('Data in PC space (top 3 PCs)');
imagesc(score(:,1:3));
xlabel('PC score (coordinate)');
ylabel('trial');

figure;
plot(explained,'-o');
xlabel('PC');
ylabel('% variance explained');


%% plot PCs
%

for condition = metadata.contextRoles

    figure;

    hold on;

    for who = metadata.subjects
        
        for run = 1:metadata.runsPerSubject
            which = which_trials & data.isTrain & data.runId == run & strcmp(data.contextRole, condition) & strcmp(data.participant, who);
            if sum(which) == 0, continue; end % only blocks in that condition
            assert(sum(which) == metadata.trainingTrialsPerRun);
            
            
            plot(score(which, 1), score(which, 2), '.-');
%            break;
        end
        
        break;
    end

    title(condition{1});
    hold off;
end

