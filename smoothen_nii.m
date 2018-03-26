function smoothen_nii(filename, how)
% Smoothen a .nii file and display it
% how = 'gauss' -> 3D Gaussian filter
% how = 'dilate' -> spherical sliding window max()
%
if ~exist('how', 'var')
    how = 'dilate';
end

[~, V, Y] = load_mask(filename);
V.fname = 'temp.nii'; % change immediately!

Y(Y < 4) = 0; % threshold

switch how
    case 'gauss'
        sigma = 1;
        scale = 3;
        B = imgaussfilt3(Y, sigma) * scale;
    case 'dilate'
        r = 3;
        B = Y;
        B(isnan(B)) = 0;
        s = strel('sphere', r);
        B = imdilate(B, s);
    otherwise
        assert(false);
end
 


spm_write_vol(V,B);
bspmview(V.fname, 'masks/mean.nii');
