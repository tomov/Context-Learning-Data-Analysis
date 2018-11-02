inds = 1:10;

N = ccnl_searchlight_rdms(context_expt(), 1, inds, getGoodSubjects());


[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());
which_rows = data.which_rows & data.isTrain; % Look at training trials only

whole_brain = load_mask('masks/mask.nii');
[x, y, z] = ind2sub(size(whole_brain), find(whole_brain)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
x = x(inds);
y = y(inds);
z = z(inds);

r = 4 / 1.5;

S = rdms_get_searchlight(data, metadata, which_rows, x, y, z, r, true, false, true); % use pregen'd betas, use tmaps, use nosmooth

n = length(getGoodSubjects());
for i = 1:length(inds)
    for s = 1:n
        a = N(i).subj(s).RDM;
        b = S(i).RDMs(:,:,s);
       
        assert(N(i).num_voxels == S(i).n_voxels);
        assert(isequal(a, b));
        %assert(immse(a(:), b(:)) < 1e-9);
    end
end
