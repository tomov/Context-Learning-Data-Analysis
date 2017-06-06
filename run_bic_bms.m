% For each mask (ROI),
%     computes BIC for each model for each subject
%     populates the log model evidences array (rows = subjects, columns =
%     models) which gets passed to BMS
%     which estimates the excedence probabilities (????)
% 
%

 % which ROIs to look at (one at a time) as paths to nifti files
roi_masks = ...
        {fullfile('masks', 'mask.nii'), ...
         fullfile('masks', 'hippocampus.nii'), ...
         fullfile('masks', 'ofc.nii'), ...
         fullfile('masks', 'striatum.nii'), ...
         fullfile('masks', 'vmpfc.nii'), ...
         fullfile('masks', 'bg.nii'), ...
         fullfile('masks', 'pallidum.nii'), ...
         fullfile('masks', 'v1.nii'), ...
         fullfile('masks', 'visual.nii'), ...
         fullfile('masks', 'motor.nii'), ...
         fullfile('masks', 'sensory.nii')};
roi_masknames = {};
for i = 1:numel(roi_masks)
    [~, roi_masknames{i}, ~] = fileparts(roi_masks{i});
end
     
% which subjects to analyze (as indices of the subjects array returned by contextGetSubjectsDirsAndRuns)
subjects = getGoodSubjects();

% which models to consider (as the glmodel value passed to context_create_multi)
%glmodels = [141 136 138]; % all trials with trial_onset regressor
glmodels = [1 140 139]; % training trials only
%models = [9 27 8];

alphas = [];
exp_rs = [];
xps = [];
pxps = [];
bors = [];
bics = [];
for roi=roi_masks
    lme = []; % log model evidence
    row_bics = [];
    for glmodel=glmodels
        bic = ccnl_bic(context_expt(), glmodel, roi{1}, subjects);
        lme = [lme, -0.5 * bic];
        row_bics = [row_bics, bic];
    end
    bics = [bics; row_bics];
        
    [alpha,exp_r,xp,pxp,bor] = bms(lme);
    alphas = [alphas; alpha];
    exp_rs = [exp_rs; exp_r];
    xps = [xps; xp];
    pxps = [pxps; pxp];
    bors = [bors; bor];
end

roi_masknames{1} = 'whole_brain';
T = array2table(xps, 'RowNames', roi_masknames, 'VariableNames', {'Full_GLM_feedback_onset', 'Full_GLM_trial_onset', 'Null_GLM'});
