% bspview -> save -> save current cluster (binary mask)
% creates a .hdr file and a .img file
% it's nice to have them coalesced into a single .nii file so we can load
% it with spm_read.
% use this script to do just that.
% also keep in mind this is super customized to remove voxels from the left
% hemisphere. you might want to comment that stuff out
%


% to convert from .hdr and .img (as bspmview outputs them) to nii
%
%prefix = 'masks/ClusterMask_spmT_0001_x=32_y=-66_z=50_1172voxels';
%prefix = 'masks/ClusterMask_spmT_0001_x=32_y=-66_z=50_458voxels';
%prefix = 'masks/ClusterMask_spmT_0001_x=32_y=-66_z=50_192voxels';
%prefix = 'masks/ClusterMask_spmT_0001_x=32_y=-66_z=50_59voxels';
%prefix = 'masks/light_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_911voxels';
%prefix = 'masks/light_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_290voxels';
%prefix = 'masks/light_ClusterMask_searchlight_tmap_x=36_y=-58_z=44_60voxels';
%prefix = 'masks/light_prior_trial-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-54_y=12_z=26_315voxels';
%prefix = 'masks/light_prior_trial-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=46_y=12_z=22_279voxels';
%prefix = 'masks/light_prior_trial-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=46_y=14_z=22_28voxels';
%prefix = 'masks/light_prior_trial-onset_L_dlPFC1_ClusterMask_searchlight_tmap_x=-36_y=10_z=24_17voxels';
%prefix = 'masks/light_prior_trial-onset_L_dlPFC2_ClusterMask_searchlight_tmap_x=-56_y=14_z=26_34voxels';

%prefix = 'masks/light_posterior_feedback-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=54_y=24_z=34_117voxels';
%prefix = 'masks/light_posterior_feedback-onset_R_dlPFC_ClusterMask_searchlight_tmap_x=54_y=24_z=34_1989voxels';
%prefix = 'masks/light_posterior_feedback-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-50_y=24_z=26_1038voxels';
%prefix = 'masks/light_posterior_feedback-onset_L_dlPFC_ClusterMask_searchlight_tmap_x=-52_y=22_z=26_61voxels';
%prefix = 'masks/light_posterior_feedback-onset_ACC_ClusterMask_searchlight_tmap_x=4_y=22_z=44_253voxels';
%prefix = 'masks/light_posterior_feedback-onset_ACC_ClusterMask_searchlight_tmap_x=6_y=22_z=44_37voxels';
%prefix = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=36_602voxels';
%prefix = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=34_119voxels';
%prefix = 'masks/light_posterior_feedback-onset_L_AG_ClusterMask_searchlight_tmap_x=-26_y=-60_z=34_10voxels';

prefix = 'masks/..';

nii = load_nii([prefix, '.hdr']);
save_nii(nii, [prefix, '.nii']);

% to trim the bad voxels in the wrong hemisphere
%
[~, V, Y] = load_mask([prefix, '.nii']);
V.fname = [prefix, '_edited.nii'];

%{
for x = 1:size(Y, 1)
    for y = 1:size(Y, 2)
        for z = 1:size(Y, 3)
            mni = cor2mni([x y z], V.mat);
            if mni(1) > 0
                Y(x,y,z) = 0;
            end
        end
    end
end
%}

spm_write_vol(V, Y);

check_mask(V.fname);

