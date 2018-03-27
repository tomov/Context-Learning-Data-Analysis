% RUN using might_group.m

%% classifier results
%

% ------ z-run ----------

% aggregate
bspmview('might/gnb_searchmight_accuracy_countmap_trial_onset_folds=3_r=2.6667_z-run_use_nosmooth=1_use_tmaps=0.nii', 'masks/mean.nii');
% aggregate, smooth w/ spheres  
smoothen_nii('might/gnb_searchmight_accuracy_countmap_trial_onset_folds=3_r=2.6667_z-run_use_nosmooth=1_use_tmaps=0.nii')

% individual subject
bspmview('might/gnb_searchmight_accuracy_trial_onset_subj=6_folds=3_r=2.6667_z-run_use_nosmooth=1_use_tmaps=0_alpha=0.050.nii', 'masks/mean.nii');


% ------- z-none ---------
% aggregate
bspmview('might/gnb_searchmight_accuracy_countmap_trial_onset_folds=3_r=2.6667_z-none_use_nosmooth=1_use_tmaps=0.nii', 'masks/mean.nii'); 
% aggregate, smooth w/ spheres -- note more stuff b/c extent = 5 in bspmview above
smoothen_nii('might/gnb_searchmight_accuracy_countmap_trial_onset_folds=3_r=2.6667_z-none_use_nosmooth=1_use_tmaps=0.nii')

% feedback onset
bspmview('might/gnb_searchmight_accuracy_countmap_feedback_onset_folds=3_r=2.6667_z-none_use_nosmooth=1_use_tmaps=0.nii', 'masks/mean.nii'); 
smoothen_nii('might/gnb_searchmight_accuracy_countmap_feedback_onset_folds=3_r=2.6667_z-none_use_nosmooth=1_use_tmaps=0.nii')


% QUESTIONS:
% - feedback onset vs. trial onset
% - visualizing classifier results
% - re-doing RSA with z-run & nosmooth?
% - z-run okay? vs. z-run-voxel?
% - correlating with behavior -- for each subject, sphere around peak accuracy voxel in each GLM cluster?
