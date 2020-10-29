

% create EXPTs for leave-one-run-out CV

test_partitions = {logical([1 1 1 0 0 0 0 0 0]), 
                   logical([0 0 0 1 1 1 0 0 0]),
                   logical([0 0 0 0 0 0 1 1 1])};

goodSubjs = getGoodSubjects();

for k = 1:3
    test = test_partitions{k};
    train = ~test_partitions{k};

    EXPT_train{k} = context_expt();
    EXPT_train{k}.modeldir = fullfile(EXPT_train{k}.exptdir, sprintf('glmOutput_CV_k=%d_train', k));
    for s = 1:length(goodSubjs)
        subj = goodSubjs(s);
        EXPT_train{k}.subject(subj).functional = EXPT_train{k}.subject(subj).functional(train);
    end

    EXPT_test{k} = context_expt();
    EXPT_test{k}.modeldir = fullfile(EXPT_test{k}.exptdir, sprintf('glmOutput_CV_k=%d_test', k));
    for s = 1:length(goodSubjs)
        subj = goodSubjs(s);
        EXPT_test{k}.subject(subj).functional = EXPT_test{k}.subject(subj).functional(test);
    end
end
