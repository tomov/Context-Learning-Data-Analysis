function masks = isl_create_masks(do_create, sphere)

    % table 3 in Jneuro 2018 paper
mni = [ 44 12 50; ...
        56 -60 30; ...
        -38 26 26; ...
        24 64 -14; ...
        -40 56 6; ...
        -40 -62 -46; ...
        2 38 48; ...
        62 -44 -12; ...
        -34 -66 -46; ...
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
