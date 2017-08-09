function [atlaslabels, atlas0] = bspm_setatlas(tmap_filename, atlas_dirpath, atlas_name)
% setatlas from bspmview
%
    atlas_vol = fullfile(atlas_dirpath, sprintf('%s_Atlas_Map.nii', atlas_name));
    if ~exist(atlas_vol, 'file')
        atlas_vol = fullfile(atlas_dirpath, sprintf('%s_Atlas_Map.nii.gz', atlas_name));
    end
    atlas_labels    = fullfile(atlas_dirpath, sprintf('%s_Atlas_Labels.mat', atlas_name));
    int             = 0; % 0=Nearest Neighbor, 1=Trilinear(default)
    atlasvol        = bspm_reslice_image(atlas_vol, tmap_filename, int);
    atlasvol        = single(round(atlasvol(:)))';
    load(atlas_labels);
    atlaslabels   = atlas;
    atlas0        = atlasvol;
