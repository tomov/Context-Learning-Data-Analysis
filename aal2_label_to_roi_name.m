function roi = aal2_label_to_roi_name(roi_label, mni)

% Function that converts AAL2 labels from e.g. bspmview tables into real anatomical ROI names.
% 
% INPUT:
% roi_label = ROI label as output by bspmview, e.g. 'Angular_R'.
% mni = optional MNI coordinates to include laterality
%
% OUTPUT:
% roi = the actual anatomical name, e.g. 'Angular gyrus' or 'Angular gyrus (R)'. If not round, returns []
%

% Anatomical region name, AAL2 label
% Table 2 from Rolls et al., Implementation of a new parcellation of the orbitofrontal cortex in the
% automated anatomical labeling atlas (NeuroImage, 2015)
%
aal2_labels = {...
'Precentral gyrus',  'Precentral' ;
'Postcentral gyrus',  'Postcentral' ;
'Rolandic operculum',  'Rolandic_Oper' ;
'Superior frontal gyrus, dorsolateral',  'Frontal_Sup' ;
'Middle frontal gyrus',  'Frontal_Mid' ;
'IFG pars opercularis',  'Frontal_Inf_Oper' ;
'IFG pars triangularis',  'Frontal_Inf_Tri' ;
'Superior frontal gyrus, medial',  'Frontal_Sup_Med' ;
'Supplementary motor area',  'Supp_Motor_Area' ;
'Paracentral lobule',  'Paracentral_Lobule' ;
'Superior frontal gyrus, medial orbital',  'Frontal_Med_Orb' ;
'IFG pars orbitalis',  'Frontal_Inf_Orb' ;
'Gyrus rectus',  'Rectus' ;
'Medial orbital gyrus',  'OFCmed' ;
'Anterior orbital gyrus',  'OFCant' ;
'Posterior orbital gyrus',  'OFCpost' ;
'Lateral orbital gyrus',  'OFClat' ;
'Olfactory cortex',  'Olfactory' ;
'Superior temporal gyrus',  'Temporal_Sup' ;
'Heschl''s gyrus',  'Heschl' ;
'Middle temporal gyrus',  'Temporal_Mid' ;
'Inferior temporal gyrus',  'Temporal_Inf' ;
'Superior parietal gyrus',  'Parietal_Sup' ;
'Inferior parietal gyrus, excluding supramarginal and angular gyri',  'Parietal_Inf' ;
'Angular gyrus',  'Angular' ;
'Supramarginal gyrus',  'SupraMarginal',;
'Precuneus',  'Precuneus' ;
'Superior occipital gyrus',  'Occipital_Sup' ;
'Middle occipital gyrus',  'Occipital_Mid' ;
'Inferior occipital gyrus',  'Occipital_Inf' ;
'Cuneus',  'Cuneus' ;
'Calcarine fissure and surrounding cortex',  'Calcarine' ;
'Lingual gyrus',  'Lingual' ;
'Fusiform gyrus',  'Fusiform' ;
'Temporal pole: superior temporal gyrus',  'Temporal_Pole_Sup' ;
'Temporal pole: middle temporal gyrus',  'Temporal_Pole_Mid' ;
'Anterior cingulate \& paracingulate gyri',  'Cingulate_Ant' ;
'Middle cingulate \& paracingulate gyri',  'Cingulate_Mid' ;
'Posterior cingulate gyrus',  'Cingulate_Post' ;
'Hippocampus',  'Hippocampus';
'Parahippocampal gyrus',  'ParaHippocampal' ;
'Insula',  'Insula' ;
'Amygdala',  'Amygdala' ;
'Caudate nucleus',  'Caudate' ;
'Lenticular nucleus, Putamen',  'Putamen' ;
'Lenticular nucleus, Pallidum',  'Pallidum' ;
'Thalamus',  'Thalamus';
'Cerebellum', 'Cerebelum'};

roi = roi_label;

for j=1:size(aal2_labels, 1)
    if startsWith(roi_label, aal2_labels{j, 2})
        roi = aal2_labels{j, 1};
        if exist('mni', 'var') % optionally include laterality
            if mni(1) < 0
                roi = [roi, ' (L)'];
            else
                roi = [roi, ' (R)'];
            end
        end
        break;
    end
end
