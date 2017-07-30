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

nii = load_nii([prefix, '.hdr']);
save_nii(nii, [prefix, '.nii']);

% to trim the bad voxels in the wrong hemisphere
%
[~, V, Y] = load_mask([prefix, '.nii']);
V.fname = [prefix, '_edited.nii'];

for x = 1:size(Y, 1)
    for y = 1:size(Y, 2)
        for z = 1:size(Y, 3)
            mni = cor2mni([x y z], V.mat);
            if mni(1) < 0
                Y(x,y,z) = 0;
            end
        end
    end
end

spm_write_vol(V, Y);

check_mask(V.fname);

