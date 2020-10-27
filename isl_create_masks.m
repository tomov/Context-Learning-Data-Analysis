function masks = isl_create_masks(do_create, sphere)

mni = [ 44 20 -6; ... % AI  context
        -46 -38 44; ... % IPL  context
        42 24 24; ... % IFG  Agency
        36 16 6; ... % AI  agency
        42 28 26; ... % IFG  vgdl
        -30 28 2; ... % AI vgdl
        -50 44 12; ... % IFG  vgdl
        ];


r = sphere / 1.5; % 10 mm sphere

[mask, V, Y] = load_mask('masks/mask.nii');
V.fname = 'masks/isl.nii'; % change immediately!

for i = 1:size(mni, 1)
    cor = mni2cor(mni(i,:), V.mat);

    V.fname = fullfile('masks', sprintf('isl_%d_%d_%d_r=%.1fmm.nii', mni(i,1), mni(i,2), mni(i,3), sphere));

    if do_create
        [sphere_mask] = ccnl_create_spherical_mask(cor(1), cor(2), cor(3), r);
        sphere_mask = sphere_mask & mask; % exclude voxels outside of brain
        spm_write_vol(V, sphere_mask);
    end

    masks{i} = V.fname;
end

%check_mask('masks/badre_rlpfc.nii');
