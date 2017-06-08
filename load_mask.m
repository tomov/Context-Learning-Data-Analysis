function mask = load_mask(mask)
% Load a .nii file as a SPM mask
%
Vmask = spm_vol(mask);
Ymask = spm_read_vols(Vmask);
mask = Ymask~=0 & ~isnan(Ymask);