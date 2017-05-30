function check_mask( mask )
    % check a mask against the mean structural, make sure it looks right
    %
    EXPT = context_expt();
    P = {fullfile(EXPT.modeldir, 'mean.nii'), fullfile(mask)};
    spm_check_registration(char(P));
end

