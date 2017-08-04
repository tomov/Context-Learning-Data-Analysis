function [C, CI, region, extent, stat, mni, cor, results_table] = bspm_extract_clusters(tmap_filename, p, direct, alpha, Dis, Num)
%
% Extract clusters and peak voxels from a t-map contrast AFTER cluster FWE
% correction. Exactly the same as bspmview -- uses the same functions. As a
% bonus, even spits out a results table. Pretty cool huh
%
% INPUT:
% tmap_filename = path to .nii file with the contrast t-map, e.g.
%                 '../neural/model154/con10/spmT_0001.nii'
% p = p-value threshold for the individual voxels, e.g. 0.001
% direct = sign of the activations to look at; should be one of +, -, or
%          +/-
% alpha = significance level for cluster FWE correction, default 0.05 in bspmview
% Dis = separation for cluster maxima, default 20 in bspmview
% Num = numpeaks for cluster maxima, default 3 in bspmview
% 
% OUTPUT:
% C = volume with cluster size for each voxel
% CI = volume with cluster index for each voxel <-- that's the name of the
%      game; 
% region = labels for the peak voxels
% extent = size of the cluster
% stat = t-statistic for the peak voxels
% mni = MNI coordinates of peak voxels
% cor = coordinates of peak voxel in native space (can be used as indices
%       in the C and CI volumes)
% results_table = what Save Results Table in bspmview would spit out 
% 

[~, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

%tmap_filename = '../neural/model154/con10/spmT_0001.nii';
%p = 0.001; % p-value threshold for individual voxels
%direct = '+/-'; % positive or negative activations
%alpha = 0.001;
%Dis = 20;
%Num = 1;


df = metadata.N - 1; % = degrees of freedom for t-test = # of subjects - 1

atlas_dirpath = '/Users/momchil/Dropbox/Research/libs/bspmview/supportfiles';
atlas_name = 'AAL2';

assert(ismember(direct, {'+', '-', '+/-'}));


% get cluster extent threshold
%
extent_thresh = bspm_cluster_correct(tmap_filename, df, direct, p, alpha);


% get cluster indices
%
thresh = spm_invTcdf(1-p, df);
V = spm_vol(tmap_filename);
Y = spm_read_vols(V);
Y(isnan(Y)) = 0;
[clustsize, clustidx]   = bspm_getclustidx(Y, thresh, extent_thresh);
di = strcmpi({'+' '-' '+/-'}, direct);
C = clustsize(di, :);
CI= clustidx(di, :);


% set some input params
%
M           = V.mat;         %-voxels to mm matrix
DIM         = V.dim';
VOX         = abs(diag(M(:,1:3)));
[x,y,z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ0         = [x(:)';y(:)';z(:)'];
RCP         = XYZ0;
RCP(4,:)    = 1;
XYZmm0       = M(1:3,:)*RCP;


% set thresh
%
idx = find(C > 0);
if find(di) == 2
    Z     = abs(Y(idx));
else
    Z     = Y(idx);
end
Nunique   = length(unique(Z));
XYZ       = XYZ0(:,idx);
XYZmm     = XYZmm0(:,idx);


% set maxima
%
switch char(direct)
    case {'+', '-'}
        LOCMAX      = bspm_getmaxima(Z, XYZ, M, Dis, Num);
        if strcmp(direct, '-')
            LOCMAX(:,2) = -LOCMAX(:,2);
        end
    otherwise
        POS         = bspm_getmaxima(Z, XYZ, M, Dis, Num);
        NEG         = bspm_getmaxima(Z*-1, XYZ, M, Dis, Num);
        if ~isempty(NEG), NEG(:,2) = NEG(:,2)*-1; end
        LOCMAX      = [POS; NEG];
end


% get labels
%
[atlaslabels, atlas] = bspm_setatlas(tmap_filename, atlas_dirpath, atlas_name);
LABELS = bspm_getregionnames(LOCMAX(:,3:5)', atlaslabels, atlas, XYZmm0);

% generate resultstable
%
results_table = [cell(size(LABELS)) LABELS num2cell(LOCMAX)];
results_table{1,1} = 'Positive';
if any(LOCMAX(:,2)<0)
    tmpidx = find(LOCMAX(:,2)<0);
    results_table{tmpidx(1),1} = 'Negative';
end

% set outputs
%
C = reshape(C, size(Y));
CI = reshape(CI, size(Y));
region = LABELS;
extent = LOCMAX(:, 1);
stat = LOCMAX(:, 2);
mni = LOCMAX(:, 3:5);
cor = mni2cor(mni, V.mat);

%% for sanity checks w/ the real bspmview
%
%{

LOCMAX_orig = LOCMAX;

load('~/Downloads/bspm_setmaxima.mat');

assert(immse(LOCMAX_orig, LOCMAX) < 1e-5);


C_orig = C;
CI_orig = CI;


load('~/Downloads/bspm_correct.mat');
disp(unique(C));
assert(isequal(C_orig, C));

[st.ol.C0, st.ol.C0IDX] = bspm_getclustidx(st.ol.Y, thresh, extent);
C = st.ol.C0(di,:);
CI= st.ol.C0IDX(di, :);
disp(unique(C));
disp(unique(CI));
assert(isequal(C_orig, C));
assert(isequal(CI_orig, CI));


load('~/Downloads/bspm_od.mat')
[clustsize, clustidx]   = bspm_getclustidx(od, thresh, extent);
C = clustsize(di, :);
CI= clustidx(di, :);
disp(unique(C));
disp(unique(CI));
assert(isequal(C_orig, C));
assert(isequal(CI_orig, CI));

% '../neural/model154/con10/spmT_0001.nii'
%}
