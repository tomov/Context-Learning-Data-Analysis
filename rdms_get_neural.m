function Neural = rdms_get_neural(data, metadata, which_rows)

% Compute the neural RDMs
% Normalized correlation for neural data for different ROIs
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
%
% OUTPUT:
% Neural = struct array of RDMs

neural_idx = 0;

disp('Computing neural RDMs...');
tic

% for each ROI, take betas at both trial onset and feedback onset
%
for event = {'trial_onset', 'feedback_onset'}
    event = event{1};

    % Load neural data
    %
    whole_brain_betas = get_betas('masks/mask.nii', event, data, metadata);

    hippocampus_mask = load_mask('masks/hippocampus.nii');
    hippocampus_betas = get_betas_submask(hippocampus_mask, whole_brain_betas);
    [hippocampusRDMs, avgHippocampusRDM] = compute_rdms(hippocampus_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = hippocampusRDMs;
    Neural(neural_idx).RDM = avgHippocampusRDM;
    Neural(neural_idx).name = ['hippocampus_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    ofc_mask = load_mask('masks/ofc.nii');
    ofc_betas = get_betas_submask(ofc_mask, whole_brain_betas);
    [ofcRDMs, avgOfcRDM] = compute_rdms(ofc_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = ofcRDMs;
    Neural(neural_idx).RDM = avgOfcRDM;
    Neural(neural_idx).name = ['OFC_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    med_ofc_mask = load_mask('masks/med_ofc.nii');
    med_ofc_betas = get_betas_submask(med_ofc_mask, whole_brain_betas);
    [medOfcRDMs, avgMedOfcRDM] = compute_rdms(med_ofc_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = medOfcRDMs;
    Neural(neural_idx).RDM = avgMedOfcRDM;
    Neural(neural_idx).name = ['mOFC_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    vmpfc_mask = load_mask('masks/vmpfc.nii');
    vmpfc_betas = get_betas_submask(vmpfc_mask, whole_brain_betas);
    [vmpfcRDMs, avgVmpfcRDM] = compute_rdms(vmpfc_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = vmpfcRDMs;
    Neural(neural_idx).RDM = avgVmpfcRDM;
    Neural(neural_idx).name = ['vmPFC_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    striatum_mask = load_mask('masks/striatum.nii');
    striatum_betas = get_betas_submask(striatum_mask, whole_brain_betas);
    [striatumRDMs, avgStriatumRDM] = compute_rdms(striatum_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = striatumRDMs;
    Neural(neural_idx).RDM = avgStriatumRDM;
    Neural(neural_idx).name = ['Striatum_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    pallidum_mask = load_mask('masks/pallidum.nii');
    pallidum_betas = get_betas_submask(pallidum_mask, whole_brain_betas);
    [pallidumRDMs, avgPallidumRDM] = compute_rdms(pallidum_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = pallidumRDMs;
    Neural(neural_idx).RDM = avgPallidumRDM;
    Neural(neural_idx).name = ['Pallidum_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    v1_mask = load_mask('masks/v1.nii');
    v1_betas = get_betas_submask(v1_mask, whole_brain_betas);
    [v1RDMs, avgV1RDM] = compute_rdms(v1_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = v1RDMs;
    Neural(neural_idx).RDM = avgV1RDM;
    Neural(neural_idx).name = ['V1_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    m1_mask = load_mask('masks/m1.nii');
    m1_betas = get_betas_submask(m1_mask, whole_brain_betas);
    [m1RDMs, avgM1RDM] = compute_rdms(m1_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = m1RDMs;
    Neural(neural_idx).RDM = avgM1RDM;
    Neural(neural_idx).name = ['M1_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    s1_mask = load_mask('masks/s1.nii');
    s1_betas = get_betas_submask(s1_mask, whole_brain_betas);
    [s1RDMs, avgS1RDM] = compute_rdms(s1_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = s1RDMs;
    Neural(neural_idx).RDM = avgS1RDM;
    Neural(neural_idx).name = ['S1_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    fusiform_mask = load_mask('masks/fusiform.nii');
    fusiform_betas = get_betas_submask(fusiform_mask, whole_brain_betas);
    [fusiformRDMs, avgFusiformRDM] = compute_rdms(fusiform_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = fusiformRDMs;
    Neural(neural_idx).RDM = avgFusiformRDM;
    Neural(neural_idx).name = ['Fusiform_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    angular_mask = load_mask('masks/angular.nii');
    angular_betas = get_betas_submask(angular_mask, whole_brain_betas);
    [angularRDMs, avgAngularRDM] = compute_rdms(angular_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = angularRDMs;
    Neural(neural_idx).RDM = avgAngularRDM;
    Neural(neural_idx).name = ['Angular_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    mid_front_mask = load_mask('masks/mid_front.nii');
    mid_front_betas = get_betas_submask(mid_front_mask, whole_brain_betas);
    [midFrontRDMs, avgMidFrontRDM] = compute_rdms(mid_front_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = midFrontRDMs;
    Neural(neural_idx).RDM = avgMidFrontRDM;
    Neural(neural_idx).name = ['MidFront_', event(1)];
    Neural(neural_idx).color = [0 1 0];

    dl_sup_front_mask = load_mask('masks/dl_sup_front.nii');
    dl_sup_front_betas = get_betas_submask(dl_sup_front_mask, whole_brain_betas);
    [dlSupFrontRDMs, avgDlSupFrontRDM] = compute_rdms(dl_sup_front_betas, 'cosine', data, metadata, which_rows);
    neural_idx = neural_idx + 1;
    Neural(neural_idx).RDMs = dlSupFrontRDMs;
    Neural(neural_idx).RDM = avgDlSupFrontRDM;
    Neural(neural_idx).name = ['dlSupFront_', event(1)];
    Neural(neural_idx).color = [0 1 0];
end

disp('Computed neural RDMs.');
toc
