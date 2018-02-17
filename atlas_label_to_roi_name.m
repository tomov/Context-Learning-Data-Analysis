function roi = atlas_label_to_roi_name(atlas_name, roi_label, mni)

% Function that converts AAL2 labels from e.g. bspmview tables into real anatomical ROI names.
% 
% INPUT:
% atlas_name = name of the atlas to use
% roi_label = ROI label from that atlas, e.g. 'Angular_R' for AAL2
% mni = MNI coordinates to infer laterality
%
% OUTPUT:
% roi = the actual anatomical name, e.g. 'Angular gyrus (R)'. If not round,
%       returns the passed roi_label
%

roi = roi_label;

switch atlas_name
    case 'AAL2'
        roi = aal2_label_to_roi_name(roi, mni);

    case 'AnatomyToolbox'
        hemi = '';
        if startsWith(roi, 'L ')
            roi = roi(3:end);
            hemi = 'L';
        elseif startsWith(roi, 'R ')
            roi = roi(3:end);
            hemi = 'R';
        end

        space = find(roi == ' ' | roi == '-');
        if ~isempty(space)
            roi = [roi(1:space), lower(roi(space+1:end))];
        end

        if ~isempty(hemi)
            roi = [roi, ' (', hemi, ')'];
        end

    case 'HarvardOxford-maxprob-thr0' 
        hemi = '';
        if  mni(1) < 0
            roi = [roi, ' (L)'];
        else
            roi = [roi, ' (R)'];
        end
        
    case 'Talairach'
        roi = roi;

    case 'Brodmann' 
        hemi = '';
        if  mni(1) < 0
            roi = [roi, ' (L)'];
        else
            roi = [roi, ' (R)'];
        end
        
    otherwise
        assert(false, 'atlas_name be one of the above');
end
