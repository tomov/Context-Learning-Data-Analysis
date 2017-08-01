function [clustsize, clustidx]   = bspm_getclustidx(rawol, u, k)

% Get cluster indices after running bspm_cluster_correct(). Shamelessly stolen from bspmview
%
% INPUT:
% im = raw volumes from the tmap .nii file, thresholded appropriately with threshold u
% u = p-value threshold that was passed to bspm_cluster_correct(), e.g. 0.001
% k = cluster extent threshold (FWE or FDR), result from bspm_cluster_correct()
%
% OUTPUT:
% 

% raw data to XYZ
DIM         = size(rawol);
[X,Y,Z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ         = [X(:)';Y(:)';Z(:)'];
pos         = zeros(1, size(XYZ, 2));
neg         = pos;
clustidx    = zeros(3, size(XYZ, 2));

% positive
supra = (rawol(:)>u)';
if sum(supra)
    tmp         = spm_clusters(XYZ(:, supra));
    clbin       = repmat(1:max(tmp), length(tmp), 1)==repmat(tmp', 1, max(tmp));
    pos(supra)  = sum(repmat(sum(clbin), size(tmp, 2), 1) .* clbin, 2)';
    clustidx(1,supra) = tmp;
end
pos(pos < k)    = 0;

% negative
rawol = rawol*-1;
supra = (rawol(:)>u)';
if sum(supra)
    tmp      = spm_clusters(XYZ(:, supra));
    clbin      = repmat(1:max(tmp), length(tmp), 1)==repmat(tmp', 1, max(tmp));
    neg(supra) = sum(repmat(sum(clbin), size(tmp, 2), 1) .* clbin, 2)';
    clustidx(2,supra) = tmp;
end
neg(neg < k) = 0;

% both
clustsize       = [pos; neg];
clustsize(3,:)  = sum(clustsize);
clustidx(3,:)   = sum(clustidx);
