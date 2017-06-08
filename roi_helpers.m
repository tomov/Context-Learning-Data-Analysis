% Helper script to generate nifti mask for given ROIs based on the AAL
% paths are relative based on Momchil's directory tree
% Requires:
%   NiFTi library: https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
%   AAL SPM toolbox: http://www.gin.cnrs.fr/AAL  (only need the AAL.nii & AAL.xlm files)
%   xml2struct script: https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
%

% Some pre-defined ROIs for your convenience
%
hippocampus_rois = {'Hippocampus_L', 'Hippocampus_R'};
ofc_rois = {'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', ...
           'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', ...
           'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R'}; 
vmpfc_rois = {'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', ...
              'Rectus_L', 'Rectus_R'};
striatum_rois = {'Caudate_L', 'Caudate_R', ...
                 'Putamen_L', 'Putamen_R'};
pallidum_rois = {'Pallidum_L', 'Pallidum_R'};
bg_rois = [striatum_rois, pallidum_rois];
occipital_rois = {'Occipital_Sup_L', 'Occipital_Sup_R', ...
                  'Occipital_Mid_L', 'Occipital_Mid_R', ...
                  'Occipital_Inf_L', 'Occipital_Inf_R'};
v1_rois = {'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', ...
               'Lingual_L', 'Lingual_R'};
visual_rois = [v1_rois, occipital_rois];
motor_rois = {'Precentral_L', 'Precentral_R'};
sensory_rois = {'Postcentral_L', 'Postcentral_R'};
rlpfc_rois = {'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R'}; % TODO wrong

interesting_rois = [v1_rois];

negative = false;

% which regions to include (must match names exactly as in AAL).
% Consult http://qnl.bu.edu/obart/explore/AAL/ or roi_names below
%
my_rois = interesting_rois;
filename = 'v1.nii'; % where to save the mask 


atlas = load_untouch_nii('../../libs/aal/atlas/AAL.nii'); % load atlas

roi_xml = xml2struct('../../libs/aal/atlas/AAL.xml'); % load annotations

% Find which values in the atlas correspond to our regions
%
roi_names = [];
roi_idxs = [];
my_idxs = [];
for i = 1:length(roi_xml.atlas.data.label)
    name = roi_xml.atlas.data.label{i}.name.Text;
    idx = str2num(roi_xml.atlas.data.label{i}.index.Text);
    if strmatch(name, my_rois)
        my_idxs = [my_idxs, idx];
    end
    roi_names = [roi_names, {name}];
    roi_idxs = [roi_idxs, idx];
end

% all_rois = containers.Map(roi_names, roi_idxs);

% Optionally flip to a negative mask
%
if ~negative
    % normal
    which = ismember(atlas.img, my_idxs);
else
    % negative -- mask == everything that's not in our ROIs
    which = ~ismember(atlas.img, my_idxs) & atlas.img > 0;
end
atlas.img(~which) = 0; % remove all regions we don't care about
atlas.img(which) = 1; % make our regions uniform

% Save the mask
%
save_untouch_nii(atlas, filename);

% Cross-register with the mean structural image to make sure the mask
% actually makes sense
%
P = {};
EXPT = context_expt();
P{1} = fullfile(filename);
P{2} = fullfile(EXPT.modeldir,'mean.nii');
spm_check_registration(char(P));