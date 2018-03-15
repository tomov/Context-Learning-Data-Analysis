% look at PCA for different size spherical ROIs; see when the Gaussian blurring stops messing things up
%


%[data,metadata,simulated] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, [1 1 0 1 0]);
%whole_brain_betas = get_betas('masks/mask.nii', 'feedback_onset', data, metadata, false);


%masks = {'masks/glm0_light_sphere_t=5.435_extent=24_roi=Frontal_Inf_Tri_L_peak=[-36_12_24].nii', ...
%         'masks/glm0_light_sphere_t=5.276_extent=27_roi=Frontal_Inf_Oper_R_peak=[48_18_26].nii', ...
%         'masks/glm0_light_sphere_t=5.687_extent=26_roi=Location not in atlas_peak=[-22_-84_-4].nii'};
%masks = masks(2);

masks = {'masks/temp.nii'};
r = 3;
create_spherical_mask(-36, 12, 24, r, 'masks/temp.nii');


betas = [];
for mask = masks
    [m, V] = load_mask(mask{1});
    b = get_activations_submask(m, whole_brain_betas);
    assert(size(b, 1) == size(data.which_rows, 1));
    betas = [betas, b];
end
assert(size(betas, 1) == size(data.which_rows, 1));


which_trials = data.which_rows; % & data.isTrain; % Look at training trials only

% 1 PCA for all subjects
%
[coeff,score,latent,tsquared,explained,mu] = pca(betas(which_trials, :));
score(which_trials, :) = score; % include dummy trials for easier indexing

% PCA for each subject
%
%clear score;
%clear coeff;
%for who_idx = 1:metadata.N
%    who = metadata.subjects{who_idx};
%    which = which_trials & strcmp(data.participant, who);
%
%    [c,s,latent,tsquared,explained,mu] = pca(betas(which, :));
%    disp(explained(1:10)');
%    score(which, :) = s; 
%    coeff{who_idx} = c;
%end


% plot PC scores over th course of the run for each subject
% rows = subjects, columns = runs
%
figure;

which_scores = 1:2;

for who_idx = 1:metadata.N
    who = metadata.subjects{who_idx};

    for run = 1:metadata.runsPerSubject

        b = score(which_trials & data.isTrain & strcmp(data.participant, who) & data.runId == run, which_scores);

        subplot(metadata.N, metadata.runsPerSubject, run + (who_idx - 1) * metadata.runsPerSubject);
        plot(b);
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        
        if who_idx == 1 && run == 1 
            legend(cellstr(num2str(which_scores')));
        end
    end
end

%{
% plot betas over th course of the run for each subject
% rows = subjects, columns = runs
%
figure;

which_betas = randsample(size(betas,2), 10); % pick 10 random betas to plot

for who_idx = 1:metadata.N
    who = metadata.subjects{who_idx};

    for run = 1:metadata.runsPerSubject

        b = betas(which_trials & data.isTrain & strcmp(data.participant, who) & data.runId == run, which_betas);

        subplot(metadata.N, metadata.runsPerSubject, run + (who_idx - 1) * metadata.runsPerSubject);
        plot(b);
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
end

%}
