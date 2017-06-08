% Script to create masks for our ROIs
%

% The ROI AAL2 labels (we make them bilateral later)
%
hippocampus = {'Hippocampus'};

ofc = {'Frontal_Med_Orb', 'Frontal_Inf_Orb_2', 'Rectus', ...
       'OFCmed', 'OFCant', 'OFCpost', 'OFClat', 'Olfactory'};
       
vmpfc = {'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', ...
              'Rectus_L', 'Rectus_R'};
          
striatum = {'Caudate_L', 'Caudate_R', ...
                 'Putamen_L', 'Putamen_R'};
             
pallidum = {'Pallidum_L', 'Pallidum_R'};

bg = [striatum, pallidum];

occipital = {'Occipital_Sup_L', 'Occipital_Sup_R', ...
                  'Occipital_Mid_L', 'Occipital_Mid_R', ...
                  'Occipital_Inf_L', 'Occipital_Inf_R'};
              
v1 = {'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', ...
               'Lingual_L', 'Lingual_R'};
           
visual = [v1, occipital];

m1 = {'Precentral_L', 'Precentral_R'};
s1 = {'Postcentral_L', 'Postcentral_R'};

% Must be equal to the variable names with the AAL2 labels
%
rois = {'hippocampus', 'ofc'};

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