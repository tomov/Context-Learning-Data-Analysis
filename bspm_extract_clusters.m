function [C, CI] = bspm_extract_clusters(tmap_filename, table_filename, threshold)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

%{
EXPT = context_expt();
mean_structural_filename = fullfile(EXPT.modeldir, 'mean.nii');
Vmean = spm_vol(mean_structural_filename);

V = spm_vol(tmap_filename);
Y = spm_read_vols(V);
V.fname = 'temp/temp.nii'; % CRUCIAL! don't overwrite t-map...

mask = Y > threshold;
spm_write_vol(V, mask);

spm_check_registration(Vmean, V);
%}

%%
%
clear all;
close all;

tmap_filename = '../neural/model154/con10/spmT_0001.nii';

p = 0.001;
alpha = 0.05;
direct = '+/-';
df = 19;
di = strcmpi({'+' '-' '+/-'}, direct);


atlas_dirpath = '/Users/momchil/Dropbox/Research/libs/bspmview/supportfiles';
atlas_name = 'AAL2';


extent = bspm_cluster_correct(tmap_filename, df, direct, p, alpha);



thresh = spm_invTcdf(1-p, df);
V = spm_vol(tmap_filename);
Y = spm_read_vols(V);
%Y(Y < thresh) = 0;
Y(isnan(Y)) = 0;

[clustsize, clustidx]   = bspm_getclustidx(Y, thresh, extent);
C = clustsize(di, :);
CI= clustidx(di, :);
disp(unique(C));
disp(unique(CI));



% read

M           = V.mat;         %-voxels to mm matrix
DIM         = V.dim';
VOX         = abs(diag(M(:,1:3)));
[x,y,z]     = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ0         = [x(:)';y(:)';z(:)'];
RCP         = XYZ0;
RCP(4,:)    = 1;
XYZmm0       = M(1:3,:)*RCP;


% set thresh

idx = find(C > 0);
if find(di) == 2
    Z     = abs(Y(idx));
else
    Z     = Y(idx);
end
Nunique   = length(unique(Z));
XYZ       = XYZ0(:,idx);
XYZmm     = XYZmm0(:,idx);
C         = C(idx);
%atlas     = st.ol.atlas0(idx);


% set maxima

%Dis = st.preferences.separation;
Dis = 20;

%Num = st.preferences.numpeaks;
Num = 3;

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
%st.ol.tab       = LOCMAX;
%st.ol.maxima    = LOCMAX(:,3:5)';

[atlaslabels, atlas] = bspm_setatlas(tmap_filename, atlas_dirpath, atlas_name);

LABELS = bspm_getregionnames(LOCMAX(:,3:5)', atlaslabels, atlas, XYZmm0);


voxels = [cell(size(LABELS)) LABELS num2cell(LOCMAX)];
voxels{1,1} = 'Positive';
if any(LOCMAX(:,2)<0)
    tmpidx = find(LOCMAX(:,2)<0);
    voxels{tmpidx(1),1} = 'Negative';
end
disp(voxels);

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
