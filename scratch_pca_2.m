% just confirming that the PC's for each condition are basically the same
%

sem = @(x) std(x) / sqrt(length(x));


[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% Load parameters
%
load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts_5_prior.mat'), 'results', 'results_options');
params = results(1).x;
options = results_options(1);

% OVERRIDE -- use pilot params as before
%
params = [0.1249 2.0064]; WRONG
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

whole_brain_betas = get_betas('masks/mask.nii', 'trial_onset', data, metadata, false);

%% pick mask
%
masks = {'masks/glm0_light_sphere_t=5.435_extent=24_roi=Frontal_Inf_Tri_L_peak=[-36_12_24].nii', ...
         'masks/glm0_light_sphere_t=5.276_extent=27_roi=Frontal_Inf_Oper_R_peak=[48_18_26].nii', ...
         'masks/glm0_light_sphere_t=5.687_extent=26_roi=Location not in atlas_peak=[-22_-84_-4].nii'};
mask = masks{2};

[m, V] = load_mask(mask);
betas = get_activations_submask(m, whole_brain_betas);
assert(size(betas, 1) == size(data.which_rows, 1));


%% PCA

which_trials = data.which_rows; % & data.isTrain; % Look at training trials only

clear coeff;
clear score;
clear latent;
clear tsquared;
clear explained;
clear mu;
for i = 1:numel(metadata.contextRoles)
    which = which_trials & strcmp(data.condition, metadata.contextRoles{i});
    [coeff{i},score{i},latent{i},tsquared{i},explained{i},mu{i}] = pca(betas(which, :));
    %score(whichs, :) = s; % include dummy trials for easier indexing
end

figure;
hold on;
for i = 1:3
    plot(coeff{i}(:,1));
end
