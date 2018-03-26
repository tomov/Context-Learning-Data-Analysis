% using the classifier posterior from various ROIs as proxy for structure posterior; see if it's better at predicting test choices than model posterior
%
% how to do it -- 
% 0) for each ROI (whole cluster) in the KL_structures GLM (maybe RSA too?)
% 1) for each subject, pick highest accuracy voxel in that cluster
% 2) use the sphere around that voxel (i.e. what got that accuracy in the first place) 
% 3) run your own classifier (to get the posterior over conditions at each time point; searchmight didn't give you that)
% 4) predict choices based on that as posterior (log lik) 
% 6) see if that's better, on average, than the model log lik




%EXPT = context_expt();
glm = 171;
contrast = 'KL_structures';

event = 'trial_onset';
r = 2.6667;

method = 'cvfitcnb';
runs = 1:9; 
trials = 11:20;
subjs = getGoodSubjects();
predict_what = 'condition';
z_score = 'z-none';

use_tmaps = false; % <-- slightly better if true; but stick with betas for consistency w/ RDMs
use_nosmooth = true; 

p = 0.001;
direct = '+';
alpha = 0.05;
Dis = 20;
Num = 1; % # peak voxels per cluster; default in bspmview is 3

if use_tmaps
    get_activations = @get_tmaps;
    load_activations = @load_tmaps;
else
    get_activations = @get_betas;
    load_activations = @load_betas;
end

%[V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_clusters(EXPT, glm, contrast, p, direct, alpha, Dis, Num);

file_format = 'might/gnb_searchmight_accuracy_%s_subj=%d_folds=3_r=%.4f_%s_use_nosmooth=%d_use_tmaps=%d.nii';
maskfile = 'masks/mask.nii';

[mask, Vmask] = load_mask(maskfile);
Vmask.fname = 'temp.nii'; % change immediately! in case we screw it up
[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
%activations = get_activations(maskfile, event, data, metadata, use_nosmooth);

%for i = 1:size(region, 1) % for each ROI
    fprintf('ROI = %s\n', region{i});

    clust_idx = CI(cor(i,1), cor(i,2), cor(i,3));
    clust_mask = CI == clust_idx;
    clust_vox = find(clust_mask);

    clust_mask = mask;
    clust_vox = find(mask);

    acc = [];
    for subj = subjs % for each subject
        filename = sprintf(file_format, event, subj, r, z_score, use_nosmooth, use_tmaps);
        [~, ~, amap] = load_mask(filename); % get accuracy map

        [~, j] = max(amap(clust_mask));
        [x,y,z] = ind2sub(size(clust_mask), clust_vox(j));
        fprintf('\nsubj = %d, max acc = %.4f\n', subj, amap(x,y,z));

        sphere_mask = create_spherical_mask_helper(mask, x, y, z, r, min_x, max_x, min_y, max_y, min_z, max_z, Vmask);

        sphere_activations = get_activations_submask(sphere_mask, activations);

        [inputs, targets, which_rows] = classify_get_inputs_and_targets_helper(runs, trials, [subj], sphere_activations, predict_what, z_score, data, metadata);
        inputs = inputs(:, sum(inputs,1) ~= 0); % clear columns of all 0's, previously NaNs b/c of use_nosmooth; we clear those out in might_getNeighbors()
        [classifier, outputs, accuracy, stats] = classify_train_helper(method, inputs, targets, runs, trials, [subj], []);

        %rid = data.runId(which_rows);
        %tid = data.trialId(which_rows);
        for run = runs
            m = mean(outputs(rid == run & tid > 10, :), 1);
            disp(m);
        end
        disp(targets(10:10:90,:));

        acc = [acc, accuracy];
    end

%end
