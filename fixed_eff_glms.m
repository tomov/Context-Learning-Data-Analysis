% run GLMs 194..197 with one set of regressors for all sessions (runs)
% SPM concatenation
% for CV

for model = 194:197

    good = getGoodSubjects();
    EXPT = context_expt();

    for s = 1:length(good)
        subj = good(s);

        ccnl_fmri_glm(context_expt(), model, subj, true);

        filename = sprintf('/Volumes/fMRI/ConLearn/glmOutput/model%d/subj%d/SPM.mat', model, subj);
        load(filename);

        modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);

        %{
        m1 = contains(SPM.xX.name, 'M1');
        m2 = contains(SPM.xX.name, 'M2');
        m3 = contains(SPM.xX.name, 'M3');

        name = {'M1', 'M2', 'M3'};
        X = [sum(SPM.xX.X(:,m1), 2) sum(SPM.xX.X(:,m2), 2) sum(SPM.xX.X(:,m3), 2)];
        %}

        %{
        figure;
        hold on;
        plot(X);
        legend({'M1', 'M2', 'M3'});
        %}

        regs = {'M1', 'M2', 'M3', 'trial_onset', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'constant'};

        X = zeros(size(SPM.xX.X,1),0);
        name = [];
        for i = 1:length(regs)
            j = contains(SPM.xX.name, regs{i});
            assert(sum(j) == 9);

            name = [name regs(i)];
            X = [X sum(SPM.xX.X(:,j), 2)];
        end

        %{
        figure;
        imagesc(SPM.xX.X);
        xticks(1:length(SPM.xX.name))
        xticklabels(SPM.xX.name);
        xtickangle(45);
        set(gca,'TickLabelInterpreter','none')

        figure;
        imagesc(X);
        xticks(1:length(name))
        xticklabels(name);
        xtickangle(45);
        set(gca,'TickLabelInterpreter','none')

        snatheu
        %}

        diff = length(SPM.xX.name) - length(name);

        SPM.xX.X = X;
        SPM.xX.name = name;
        SPM.xX.iC = 1:length(name);
        SPM.xX.iB = SPM.xX.iB - diff;

        save(filename, 'SPM');

        cdir = pwd;
        cd(modeldir);

        spm_spm(SPM);

        cd(cdir);
    end

end
