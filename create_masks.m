% Script to create masks for our ROIs
%

% The ROI AAL2 labels (we make them bilateral later)
%
hippocampus = {'Hippocampus'};
ofc = {'OFCmed', 'OFCant', 'OFCpost', 'OFClat'};
med_ofc = {'OFCmed'};
rectus = {'Rectus'};
vmpfc = {'Rectus', 'Frontal_Sup_Medial', 'Frontal_Med_Orb', 'Cingulate_Ant', 'OFCmed'};
striatum = {'Caudate', 'Putamen'};
pallidum = {'Pallidum'};
bg = [striatum, pallidum];
v1 = {'Calcarine'};
m1 = {'Precentral'};
s1 = {'Postcentral'};
fusiform = {'Fusiform'};
angular = {'Angular'};
mid_front = {'Frontal_Mid_2'};
dl_sup_front = {'Frontal_Sup_2'};

% Must be equal to the variable names with the AAL2 labels
%
rois = {'hippocampus', 'ofc', 'med_ofc', 'rectus', 'vmpfc', 'striatum', 'pallidum', 'bg', ...
        'v1', 'm1', 's1', 'fusiform', 'angular', 'mid_front', 'dl_sup_front'};

% Create the masks
%
for roi = rois
    roi = roi{1};
    
    % Get the AAL2 labels and create separate entries for left and right
    % hemisphere (i.e. bilateralize ROIs)
    %
    eval(['labels = ', roi, ';']);
    labels_L = cellfun(@(x) [x, '_L'], labels, 'UniformOutput', false);
    labels_R = cellfun(@(x) [x, '_R'], labels, 'UniformOutput', false);
    labels = [labels_L, labels_R];
    
    % Create and save the mask
    %
    create_mask(labels, fullfile('masks', [roi, '.nii']), true);
    create_mask(labels, fullfile('masks', [roi, '_unnormalized.nii']), false);
end