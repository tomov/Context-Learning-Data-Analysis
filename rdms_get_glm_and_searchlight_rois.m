function Neural = rdms_get_glm_and_searchlight_rois(data, metadata, which_rows)
% Compute the neural RDMs for a bunch of ROIs based on clusters from the GLMs
% and the searchlight analysis
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
%
% OUTPUT:
% Neural = struct array of RDMs

events = {'trial_onset', 'feedback_onset'};
use_tmaps = false;
use_nosmooth = false;

mask_idx = 0;

% GLM 123 (KL_structures) clusters
%
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_1172voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-1172';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_458voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-458';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_192voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-192';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_AG_ClusterMask_spmT_0001_x=32_y=-66_z=50_59voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-AG-49';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_dlPFC_ClusterMask_spmT_0001_x=48_y=24_z=32_1160voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-dlPFC_1160';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_dlPFC_ClusterMask_spmT_0001_x=46_y=26_z=26_409voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-dlPFC_409';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_dlPFC_ClusterMask_spmT_0001_x=46_y=26_z=28_166voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-dlPFC_166';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_dlPFC_ClusterMask_spmT_0001_x=46_y=26_z=28_27voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-dlPFC_27';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_L_dlPFC_ClusterMask_spmT_0001_x=-48_y=18_z=24_621voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-L-dlPFC_621';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_L_dlPFC_ClusterMask_spmT_0001_x=-38_y=14_z=26_73voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-L-dlPFC_73';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_rlPFC_ClusterMask_spmT_0001_x=36_y=56_z=-4_928voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-rlPFC_928';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_rlPFC_ClusterMask_spmT_0001_x=36_y=56_z=-2_90voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-rlPFC_90';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_R_rlPFC_ClusterMask_spmT_0001_x=38_y=54_z=-2_22voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-R-rlPFC_22';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_L_rlPFC_ClusterMask_spmT_0001_x=-40_y=56_z=4_1121voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-L-rlPFC_1121';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_L_rlPFC_ClusterMask_spmT_0001_x=-40_y=56_z=4_302voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-L-rlPFC_302';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_L_rlPFC_ClusterMask_spmT_0001_x=-40_y=56_z=2_115voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-L-rlPFC_115';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm123_L_rlPFC_ClusterMask_spmT_0001_x=-40_y=56_z=2_19voxels_edited.nii';
masks(mask_idx).rdm_name = 'GLM1-L-rlPFC_19';


% searchlight -- posterior @ feedback_onset clusters
%
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_AG_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_911voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-AG-911';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_AG_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_290voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-AG-290';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_AG_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_60voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-AG-60';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=54_y=24_z=34_1989voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-dlPFC-1989';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=54_y=24_z=34_117voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-R-dlPFC-117';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-50_y=24_z=26_1038voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-dlPFC-1038';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-52_y=22_z=26_61voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-dlPFC-61';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_ACC_ClusterMask_searchlight_tmap_x=4_y=22_z=44_253voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-ACC-253';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_ACC_ClusterMask_searchlight_tmap_x=6_y=22_z=44_37voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-ACC-37';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=36_602voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-AG-602';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=34_119voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-AG-119';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=34_10voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-post-f-L-AG-10';


% searchlight -- prior @ trial_onset clusters
%
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-54_y=12_z=26_315voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC-315';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=46_y=14_z=22_28voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC-28';

mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=46_y=12_z=22_279voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC-279';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_L_dlPFC1_ClusterMask_searchlight_tmap_x=-36_y=10_z=24_17voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC1-17';
mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/light_prior_trial-onset_L_dlPFC2_ClusterMask_searchlight_tmap_x=-56_y=14_z=26_34voxels_edited.nii';
masks(mask_idx).rdm_name = 'light-pri-t-L-dlPFC2-34';


mask_idx = mask_idx + 1;
masks(mask_idx).filename = 'masks/glm0_light_cluster_t=5.435_extent=17_roi=Frontal_Inf_Tri_L_peak=[-36_12_24].nii';
masks(mask_idx).rdm_name = 'FUCK ME';


Neural = rdms_get_neural(masks, events, data, metadata, which_rows, use_tmaps, use_nosmooth);
