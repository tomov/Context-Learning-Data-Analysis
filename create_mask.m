function [my_mask, my_vol] = create_mask(my_roi_labels, filename, normalize)

% Create a mask from a list of ROI labels.
%
% INPUT:
% my_roi_labels = cell array with AAL2 labels of the ROIs to include, e.g.
%                 {'Hippocampus_L', 'Hippocampus_R'}. Consult
%                 Edmund T. Rolls, Marc Joliot, Nathalie Tzourio-Mazoyer, 
%                 "Implementation of a new parcellation of the orbitofrontal cortex in the automated anatomical labeling atlas"
%                 Neuroimage, 2015
% filename = output filename where to save the .nii file
% normalize = true if you want to put the mask in the coordinate space of
%             our subject-group level mask and to exclude any voxels that
%             are not included in the subject group-level mask
%             false if you want to leave the mask in AAL2 space
%
% OUTPUT:
% my_mask = the actual mask as a 3D binary vector
% my_vol = the corresponding SPM volume variable
%

aal2_atlas_filename = fullfile('masks', 'AAL2.nii');
aal2_labels_filename =  fullfile('masks', 'AAL2.xml');
group_mask_filename = fullfile('masks', 'mask.nii'); % the subject group-level mask

% Load the AAL2 anatomical labels (indices) for each voxel
%
[~, aal2_vol, aal2_mask] = load_mask(aal2_atlas_filename);

% Load the AAL2 mapping from label index to label name
%
aal2_labels = xml2struct(aal2_labels_filename);
x = {};
for i = 1:numel(aal2_labels.atlas.data.label)
    label = aal2_labels.atlas.data.label{i}.name.Text;
    x = [x, {label}];
end
assert(numel(x) == numel(aal2_labels.atlas.data.label));
valid_labels = ismember(my_roi_labels, x);
if sum(~valid_labels) ~= 0
    disp('Some of your ROI labels are nonexistent!');
    disp(my_roi_labels(~valid_labels));
    assert(false);
end

% Load the subject group-level mask
%
[~, group_vol, group_mask] = load_mask(group_mask_filename);

% Find which label indices in the atlas correspond to our regions
%
my_roi_idxs = [];
for i = 1:length(aal2_labels.atlas.data.label)
    label = aal2_labels.atlas.data.label{i}.name.Text;
    idx = str2num(aal2_labels.atlas.data.label{i}.index.Text);
    if strmatch(label, my_roi_labels)
        my_roi_idxs = [my_roi_idxs, idx];
    end
end

% Construct the mask
%
my_vol = aal2_vol;
my_vol.fname = filename; % CRUCIAL! o/w overwrite AAL2.nii
my_mask = ismember(aal2_mask, my_roi_idxs);

% Optionally put the mask in our subject coordinate space,
% i.e. the coordinate space of the subject group-level mask / the betas
% and AND with our subjet group-level mask
%
if normalize
    [x, y, z] = ind2sub(size(my_mask), find(my_mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
    cor = mni2cor(cor2mni([x y z], my_vol.mat), group_vol.mat); % voxel coords in AAL2 space --> voxel coords in MNI space --> voxel coords in our space
    ind = sub2ind(size(group_mask), cor(:,1), cor(:,2), cor(:,3)); % voxel coords in our space --> voxel indices
    
    % Reconstruct mask in our space
    %
    my_vol = group_vol;
    my_vol.fname = filename; % CRUCIAL! o/w overwrite mask.nii
    my_mask = zeros(size(group_mask));
    my_mask(ind) = 1; % voxel indices --> binary mask
    
    % Only include voxels that are part of the subject group-level mask
    % i.e. that have non-NaN betas for all subjects
    %
    my_mask = my_mask & group_mask;
end

% Save the mask
%
spm_write_vol(my_vol, my_mask);