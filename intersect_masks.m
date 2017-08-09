function intersect_masks(filenames, outfile)
% Given a bunch of masks, intersect them
%
% INPUT:
% filenames = cell array of paths to .nii files with the masks
%
% OUTPUT:
% outfile = path to the .nii file where to save the resulting intersection
%

[outmask, V, ~] = load_mask(filenames{1});
V.fname = outfile; % change immediately!!!!

for i = 2:numel(filenames)
    filename = filenames{i};
    mask = load_mask(filename);
    outmask = outmask & mask;
end

spm_write_vol(V, outmask);
