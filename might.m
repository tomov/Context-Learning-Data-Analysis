addpath('/Users/momchil/Dropbox/Research/libs/SearchmightToolbox.Darwin_i386.0.2.5');

maskfile = 'masks/mask.nii';
event = 'feedback_onset';
use_nosmooth = true;

[mask] = load_mask(maskfile);
% load instead
%
%[meta] = createMetaFromMask(mask);
load('might/meta.mat');

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% fix neighbors according to spherical searchlight
%
[meta.voxelsToNeighbours, meta.numberOfNeighbours] = might_computeNeighborsSphere(mask, meta.colToCoord, r);

% this is unnecessary -- examples == betas
%
%{
betas = get_betas('masks/mask.nii', 'feedback_onset', data, metadata, false);
[x, y, z] = ind2sub(size(mask), find(mask)); % NOTE this works ONLY b/c of the way we got the betas in ccnl_get_betas, i.e. Y(mask) == Y(find(mask))
dimx = size(mask, 1);
dimy = size(mask, 2); 
dimz = size(mask, 3); 
dimt = size(betas, 1);
data = nan(dimx, dimy, dimz, dimt);
for i = 1:size(betas, 2)
    data(x(i),y(i),z(i),:) = betas(:,i);
end

[dimx,dimy,dimz,dimt] = size(data);
nvoxels = length(meta.indicesIn3D);

examples = zeros(dimt,nvoxels);

for t = 1:dimt

    volume = data(:,:,:,t);
    examples(t,:) = volume(meta.indicesIn3D);

end
%}

betas = get_betas(maskfile, event, data, metadata, use_nosmooth);
dimx = size(mask, 1);
dimy = size(mask, 2); 
dimz = size(mask, 3); 
nvoxels = size(betas, 2);


%examples = f


