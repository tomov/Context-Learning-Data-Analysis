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
close all;
clear all;

tmap_filename = '../neural/model154/con10/spmT_0001.nii';

p = 0.001;
alpha = 0.05;
direct = '+/-';
df = 19;
di = strcmpi({'+' '-' '+/-'}, direct);

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


load('~/Downloads/bspm_correct.mat');
disp(unique(C));
%Y = st.ol.Y;
[st.ol.C0, st.ol.C0IDX] = bspm_getclustidx(st.ol.Y, thresh, extent);
C = st.ol.C0(di,:);
disp(unique(C));


load('~/Downloads/bspm_od.mat')
[clustsize, clustidx]   = bspm_getclustidx(od, thresh, extent);
C = clustsize(di, :);
CI= clustidx(di, :);
disp(unique(C));

% '../neural/model154/con10/spmT_0001.nii'
