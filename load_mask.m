function [mask, Vmask, Ymask] = load_mask(mask_filename)
% Load a .nii file as a SPM mask
%
Vmask = spm_vol(mask_filename);
Ymask = spm_read_vols(Vmask);
mask = Ymask~=0 & ~isnan(Ymask);
