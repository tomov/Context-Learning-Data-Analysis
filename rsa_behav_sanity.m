[B, c1] = ccnl_behavioral_rdms(context_expt(), 1, getGoodSubjects());

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());
which_rows = data.which_rows & data.isTrain; % Look at training trials only

[M, c2, params, which_structures] = rdms_get_model_3(data, metadata, which_rows);

n = length(getGoodSubjects());
for s = 1:n
    a = B(1).subj(s).RDM;
    b = M(1).RDMs(:,:,s);
    
    assert(isequal(a, b));
    %assert(immse(a(:), b(:)) < 1e-9);
end
