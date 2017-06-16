

%{

method = 'cvglmnet';
mask = fullfile('masks', 'hippocampus.nii');
z_score = 'z-none';
predict_what = 'condition';
runs = [1:9];
trials = [1:24];
subjs = getGoodSubjects();

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, subjs);

leftout_run = 9;
runs(runs == leftout_run) = [];
which_rows = data.which_rows & ismember(data.trialId, trials) & data.runId ~= leftout_run;

m = containers.Map(metadata.subjects, subjs);
subj_id = cellfun(@(x) m(x), data.participant(which_rows));

foldid = data.runId(which_rows) + subj_id * (metadata.runsPerSubject + 1);
assert(numel(unique(foldid)) == numel(subjs) * (metadata.runsPerSubject - 1));
unique_foldid = unique(foldid);
foldid = arrayfun(@(x) find(unique_foldid == x), foldid);
assert(numel(unique(foldid)) == numel(subjs) * (metadata.runsPerSubject - 1));
assert(max(foldid) == numel(subjs) * (metadata.runsPerSubject - 1));


[classifier, inputs, targets, outputs, which_rows] = classify_train(method, runs, trials, subjs, mask, predict_what, z_score, foldid);

[test_inputs, test_targets, test_outputs, test_which_rows] = classify_test(method, classifier, leftout_run, trials, subjs, mask, predict_what, z_score);

%}

%{
EXPT = context_expt();
subj = 1;
model = 139;
modeldir = fullfile(EXPT.modeldir,['model',num2str(model)],['subj',num2str(subj)]);
load(fullfile(modeldir,'SPM.mat'));

[N,K] = size(SPM.xX.X)
%}


%{
o = all_outputs(all_which_rows, :);
t = all_targets(all_which_rows, :);
acc = classify_get_accuracy(ones(size(o)) / 3, t)
%}


%{
handle = figure;
        %set(handle, 'Position', [500, 500, 450, 200])
        
        subplot(2, 1, 1);

        P_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.isTrain & data.trialId == 20 & strcmp(data.contextRole, condition);
            
            P = simulated.P(which_rows, which_structures);             
            P_means = [P_means; mean(P, 1)];
        end

        bar(P_means);
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Posterior probability');
        legend({'M1', 'M2', 'M3'}, 'Position', [0.2 0.3 1 1]);
        title('Model');
    
        %
        % Figure 3B: Choice probabilities on test trials for model vs. humans
        %
        
        subplot(2, 1, 2);
        
        % TODO SAM superimpose human dots w/ error bars; rm error bars from model
        %
        % TODO SAM also SEMs -- within-subject errors: wse.m (in the dropbox)
        % TODO SAM also mytitle.m -- easier to plot left-oriented panel titles
        
        % Choice probabilities for model
        %
        human_means = [];
        human_sems = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & ~data.isTrain & strcmp(data.contextRole, condition);
            
            x1c1 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 0);
            x1c3 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 2);
            x3c1 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 0);
            x3c3 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 2);

            human_means = [human_means; mean(x1c1) mean(x1c3) mean(x3c1) mean(x3c3)];
            human_sems = [human_sems; sem(x1c1) sem(x1c3) sem(x3c1) sem(x3c3)];
        end

        hold on;
        %errorbar(human_means(1,:), human_sems(1,:));
        hold off;
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Choice probability');
        legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'}, 'Position', [0.07 0.2 1 1]);
        title('Human subjects');
%}
