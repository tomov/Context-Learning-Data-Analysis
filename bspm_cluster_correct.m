function [fwek, fdrk] = bspm_cluster_correct(im,df,direct,u,alpha)
% Cluster correction function copy-pasted from bspmview, with some modifications to 
% make it work standalone.
%
% INPUT:
% im = Path to .nii file with t-map, e.g. '../neural/model154/con10/spmT_0001.nii'
% df = degrees of freedom = # of subjects - 1 (for group-level contrast)
% direct = direction of t-values, should be one of {'+' '-' '+/-'}
% u = p-value threshold (exclude all voxels with greater p-values), default 0.001
% alpha = significance level of the cluster correction (default 0.05)
%
% OUTPUT:
% fwek = Cluster FWE extent threshold
% fdrk = Cluster FDR extent threshold
% 
% THIS IS A MODIFICATION OF A FUNCTION BY BSPMVIEW. THIS IS THE ORIGINAL DOCUMENTATION
%
% CLUSTER_CORRECT Computer extent for cluster-level correction
% USAGE: [k info] = cluster_correct(im,u,alpha,range)
%
%
% THIS IS A MODIFICATION OF A FUNCTION BY DRS. THOMAS NICHOLS AND MARKO
% WILKE, CorrClusTh.m. ORIGINAL DOCUMENTATION PASTED BELOW:
%
% Find the corrected cluster size threshold for a given alpha
% function [k,Pc] =CorrClusTh(SPM,u,alpha,guess)
% SPM   - SPM data structure
% u     - Cluster defining threshold
%         If less than zero, u is taken to be uncorrected P-value
% alpha - FWE-corrected level (defaults to 0.05)
% guess - Set to NaN to use a Newton-Rhapson search (default)
%         Or provide a explicit list (e.g. 1:1000) of cluster sizes to
%         search over.
%         If guess is a (non-NaN) scalar nothing happens, except the the
%         corrected P-value of guess is printed.
%
% Finds the corrected cluster size (spatial extent) threshold for a given
% cluster defining threshold u and FWE-corrected level alpha.
%
%_________________________________________________________________________
% $Id: CorrClusTh.m,v 1.12 2008/06/10 19:03:13 nichols Exp $ Thomas Nichols, Marko Wilke
if nargin < 4, u = .001; end
if nargin < 5, alpha = .05; end
if iscell(im), im = char(im); end

assert(ismember(direct, {'+' '-' '+/-'}));

V = spm_vol(im);
Ymask = spm_read_vols(V);
DIM = V.dim';
[X,Y,Z] = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ0 = [X(:)';Y(:)';Z(:)'];

%% Get Design Variable %%
if exist([fileparts(im) filesep 'I.mat'],'file')
    matfile = [fileparts(im) filesep 'I.mat'];
    flexflag = 1;
elseif exist([fileparts(im) filesep 'SPM.mat'],'file')
    matfile = [fileparts(im) filesep 'SPM.mat'];
    flexflag = 0;
else
    disp('Could not find an SPM.mat or I.mat variable, exiting.'); extent = []; info = []; return
end

%% Load and Compute Params %%
if flexflag % GLMFLEX

    II = load(matfile);
    try
        maskhdr = spm_vol(fullfile(fileparts(im), 'mask.nii'));
    catch
        maskhdr = spm_vol(fullfile(fileparts(im), 'mask.nii.gz'));
    end

    FWHM    = II.I.FWHM{1};
    R       = spm_resels_vol(maskhdr,FWHM)';
    VOX2RES = 1/prod(FWHM(~isinf(FWHM)));% voxels to resels

else % SPM

    load(matfile);
    R    = SPM.xVol.R;
    VOX2RES  = 1/prod(SPM.xVol.FWHM(~isinf(SPM.xVol.FWHM))); %-voxels to resels

end

%% SPM METHOD (spm_uc_clusterFDR)
thresh         = spm_invTcdf(1-u, df);
di             = find(strcmpi({'+' '-' '+/-'}, direct));
idx            = Ymask(:)~=0;
Z              = Ymask(idx)';
if di==2, Z = Z*-1; end
XYZ            = XYZ0(:,idx);
[fdrk,Pc,fwek] = spm_uc_clusterFDR(alpha, [1 df], 'T', R, 1, Z, XYZ, VOX2RES, thresh);
