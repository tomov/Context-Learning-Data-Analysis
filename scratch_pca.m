% here's why the PCA stuff is not legit
% the effect is driven by a few outliers
% and the distribution of PC1 is not Gaussian (clearly b/c of the outliers)
% -> the ANOVAs and t-tests get fooled; if you remove the outliers, the
% conditions are indistinguishable;
% a classifier is not able to distinguish on a per-run basis better than
% chance (except for the few outliers)
%

% PCA stuff based on https://github.com/tomov/neurolab/tree/master/exercise

% see the tmap from the prior RDMs (i.e. where the prior might be stored)
%bspmview('rdms/betas_smooth/searchlight_tmap_prior_trial_onset.nii', 'masks/mean.nii');

% get spherical masks from the peak voxels in the clusters in that tmap
%create_sphere_masks_from_contrast('rdms/betas_smooth/searchlight_tmap_prior_trial_onset.nii', 0, 'light', 0.001, '+', 0.05, 20, 1, 1.814);
% -> 
% 'masks/glm0_light_sphere_t=5.687_extent=26_roi=Location not in atlas_peak=[-22_-84_-4].nii'
% 'masks/glm0_light_sphere_t=5.435_extent=24_roi=Frontal_Inf_Tri_L_peak=[-36_12_24].nii'
% 'masks/glm0_light_sphere_t=5.276_extent=27_roi=Frontal_Inf_Oper_R_peak=[48_18_26].nii'
%
sem = @(x) std(x) / sqrt(length(x));

%% load behavioral stuff

% Load data
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

% Load parameters
%
load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts_5_prior.mat'), 'results', 'results_options');
params = results(1).x;
options = results_options(1);

% OVERRIDE -- use pilot params as before
%
params = [0.1249 2.0064]; WRONG
options.isFmriData = false;
options.fixedEffects = true;

disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);
% safeguards
%assert(options.isFmriData == true);
%assert(~options.fixedEffects);
assert(isequal(options.which_structures, [1 1 1 0]));
which_structures = logical([1 1 1 0]);

% Run the model with the parameters
%
simulated = simulate_subjects(data, metadata, params, which_structures);

%% load betas -- SLOW
%

whole_brain_betas = get_betas('masks/mask.nii', 'trial_onset', data, metadata, false);

%% load mask & corresponding betas
%
%bspmview('masks/prior_left_IFG.nii', 'masks/mean.nii');

%{
%[m, V] = load_mask('masks/prior_left_IFG.nii'); % left IFG for now TODO right one too
%mask = '';
mask = 'masks/glm0_light_sphere_t=5.276_extent=27_roi=Frontal_Inf_Oper_R_peak=[48_18_26].nii';
%mask = 'masks/glm0_light_sphere_t=5.687_extent=26_roi=Location not in atlas_peak=[-22_-84_-4].nii';

[m, V] = load_mask(mask);
betas = get_activations_submask(m, whole_brain_betas);
assert(size(betas, 1) == size(data.which_rows, 1));
%}
masks = {'masks/glm0_light_sphere_t=5.435_extent=24_roi=Frontal_Inf_Tri_L_peak=[-36_12_24].nii', ...
         'masks/glm0_light_sphere_t=5.276_extent=27_roi=Frontal_Inf_Oper_R_peak=[48_18_26].nii', ...
         'masks/glm0_light_sphere_t=5.687_extent=26_roi=Location not in atlas_peak=[-22_-84_-4].nii'};
masks = masks(2);
     
betas = [];
for mask = masks
    [m, V] = load_mask(mask{1});
    b = get_activations_submask(m, whole_brain_betas);
    assert(size(b, 1) == size(data.which_rows, 1));
    betas = [betas, b];
end
assert(size(betas, 1) == size(data.which_rows, 1));

% restrict to relevant trials
%
which_trials = data.which_rows; % & data.isTrain; % Look at training trials only

%{
figure;
title('Data in voxel space');
imagesc(betas(which_trials, :));
xlabel('voxel');
ylabel('trial');
%}

% run PCA for all subjects together
%
[coeff,score,latent,tsquared,explained,mu] = pca(betas(which_trials, :));
score(which_trials, :) = score; % include dummy trials for easier indexing

%{
figure;
title('Data in PC space (top 3 PCs)');
imagesc(score(:,1:3));
xlabel('PC score (coordinate)');
ylabel('trial');

figure;
plot(explained,'-o');
xlabel('PC');
ylabel('% variance explained');
%}


% PC1 over time for each condition, collapsed across subjects & blocks
%
figure;

subplot(1,3,1);
plot(explained,'-o');
xlabel('PC');
ylabel('% variance explained');

% training trials
subplot(1,3,2);
hold on;

cond_idx = 0;
for condition = metadata.contextRoles
    cond_idx = cond_idx + 1;

    pc1 = [];
    for t = 1:metadata.trainingTrialsPerRun
        which = which_trials & data.isTrain & data.trialId == t & strcmp(data.contextRole, condition);
        
        s = nan(metadata.runsPerContext, metadata.N); % average within subject
        for who_idx = 1:metadata.N
            s(:, who_idx) = score(which & strcmp(data.participant, metadata.subjects{who_idx}), 1);
        end
        assert(isequal(size(s), [metadata.runsPerContext, metadata.N]));
        %pc1 = [pc1, score(which, 1)];  WRONG -- do NOT take SEMs across
        %blocks and subjects ... instead, average first within subject
        %(i.e. across blocks), then report SEM across subjects:
        pc1 = [pc1, mean(s)'];
    end
    
    errorbar(mean(pc1), sem(pc1)); % average across subjects
end

legend(metadata.contextRoles);
hold off;
title(mask, 'interpreter', 'none');
xlabel('training trial #');
ylabel('PC 1');


%{
% test trials -- not pretty; not enough power
subplot(1,5,3);
hold on;

cond_idx = 0;
for condition = metadata.contextRoles
    cond_idx = cond_idx + 1;

    pc1 = [];
    for t = 1:metadata.testTrialsPerRun
        which = which_trials & data.isTest & data.trialId == t & strcmp(data.contextRole, condition);
        pc1 = [pc1, score(which, 1)];       
    end
    
    errorbar(mean(pc1), sem(pc1));
end

legend(metadata.contextRoles);
hold off;
title(mask, 'interpreter', 'none');
xlabel('test trial #');
ylabel('PC 1');
%}


% collapsed across time; training only
subplot(1,3,3);

pc1 = [];
for condition = metadata.contextRoles
    which = which_trials & data.isTrain & strcmp(data.contextRole, condition);
    
    s = nan(metadata.runsPerContext * metadata.trainingTrialsPerRun, metadata.N); % average within subject
    for who_idx = 1:metadata.N
        s(:, who_idx) = score(which & strcmp(data.participant, metadata.subjects{who_idx}), 1);
    end
    assert(isequal(size(s), [metadata.runsPerContext * metadata.trainingTrialsPerRun, metadata.N]));
    pc1 = [pc1, mean(s)'];
end
h = bar(mean(pc1), 'FaceColor',[0 .5 .5]);
hold on;
errorbar(mean(pc1), sem(pc1), '.', 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
hold off;
xticklabels(metadata.contextRoles);
ylabel('PC 1');



%{
% collapsed across time; test only -- not pretty; not enough power
subplot(1,5,5);

pc1 = [];
for condition = metadata.contextRoles
    which = which_trials & data.isTest & strcmp(data.contextRole, condition);
    pc1  = [pc1, score(which, 1)];
end
h = bar(mean(pc1), 'FaceColor',[0 .5 .5]);
hold on;
errorbar(mean(pc1), sem(pc1), '.', 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
hold off;
xticklabels(metadata.contextRoles);
ylabel('PC 1');
%}


%% second half vs. first half of block
%
figure;

pc1_first_half = [];
pc1_second_half = [];

pc1_first_half_notcollapsed = [];
pc1_second_half_notcollapsed = [];
for condition = metadata.contextRoles
    which = which_trials & data.isTrain & strcmp(data.contextRole, condition) & data.trialId <= 10;
    pc1_first_half_notcollapsed  = [pc1_first_half_notcollapsed, score(which, 1)];
    s = nan(metadata.runsPerContext * metadata.trainingTrialsPerRun / 2, metadata.N); % average within subject
    for who_idx = 1:metadata.N
        s(:, who_idx) = score(which & strcmp(data.participant, metadata.subjects{who_idx}), 1);
    end
    assert(isequal(size(s), [metadata.runsPerContext * metadata.trainingTrialsPerRun / 2, metadata.N]));    
    pc1_first_half  = [pc1_first_half, mean(s)'];
    
    which = which_trials & data.isTrain & strcmp(data.contextRole, condition) & data.trialId > 10;
    pc1_second_half_notcollapsed  = [pc1_second_half_notcollapsed, score(which, 1)];
    s = nan(metadata.runsPerContext * metadata.trainingTrialsPerRun / 2, metadata.N); % average within subject
    for who_idx = 1:metadata.N
        s(:, who_idx) = score(which & strcmp(data.participant, metadata.subjects{who_idx}), 1);
    end
    assert(isequal(size(s), [metadata.runsPerContext * metadata.trainingTrialsPerRun / 2, metadata.N]));    
    pc1_second_half  = [pc1_second_half, mean(s)'];
end
pc1_means = [mean(pc1_first_half); mean(pc1_second_half)];
pc1_sems = [sem(pc1_first_half); sem(pc1_second_half)];

h = bar(pc1_means);
xs = sort([h(1).XData + h(1).XOffset, h(2).XData + h(2).XOffset, h(3).XData + h(3).XOffset]);
pc1_means = pc1_means'; pc1_means = pc1_means(:);
pc1_sems = pc1_sems'; pc1_sems = pc1_sems(:);

hold on;
errorbar(xs, pc1_means(:), pc1_sems(:), '.', 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
hold off;
%xticklabels(metadata.contextRoles);
ylabel('PC 1');
legend(metadata.contextRoles);
xticklabels({'1st half', '2nd half'});


%% two-way anova -- is PC1 changing over time in the different conditions?
% rows = 1st half vs. 2nd half;
% columns = irrelevant, modulatory, addditive
pc1 = [pc1_first_half; pc1_second_half];
[p,tbl] = anova2(pc1, metadata.N);


%% one-way anova -- is PC different across conditions?
% columns = irrelevant, modulatory, addditive
pc1 = [pc1_first_half; pc1_second_half];
pc1 = [pc1_first_half_notcollapsed; pc1_second_half_notcollapsed];
[p,tbl] = anova1(pc1);



%% one-way anova -- is PC1 different across conditions in the 2nd half?
% columns = irrelevant, modulatory, addditive
[p,tbl] = anova1(pc1_second_half);




%% should we z-score? plot all the data! rows = subjects, columns = runs
% --> of course they're super correlated ... b/c we're correlating them... z o m f g 
% also, you're an idiot if you z-score for each run, then the average
% activation in each run would be 0! how would you compare across runs of
% different conditions?
%

figure;

for who_idx = 1:metadata.N
    who = metadata.subjects{who_idx};

    for run = 1:metadata.runsPerSubject

        b = betas(which_trials & data.isTrain & strcmp(data.participant, who) & data.runId == run, :);

        subplot(metadata.N, metadata.runsPerSubject, run + (who_idx - 1) * metadata.runsPerSubject);
        plot(b);
        set(gca, 'xtick', []);
        %set(gca, 'ytick', []);
    end
end






%% PC1 over time for each condition for each subject, collapsed across
% blocks
%
figure;

for who_idx = 1:metadata.N
    who = metadata.subjects{who_idx};

    subplot(4,5,who_idx);
    
    hold on;

    cond_idx = 0;
    for condition = metadata.contextRoles
        cond_idx = cond_idx + 1;

        pc1 = [];
        for t = 1:metadata.trainingTrialsPerRun
            which = which_trials & data.isTrain & data.trialId == t & strcmp(data.contextRole, condition) & strcmp(data.participant, who);
            pc1 = [pc1, score(which, 1)];
        end
        %errorbar(mean(pc1), sem(pc1)); % average across blocks
        plot(mean(pc1));
    end
    
    hold off;
    set(gca, 'xtick', []);
end




%% classifier: multi-class SVM
%

which = which_trials & data.isTrain;

% average the betas
X = [];
Y = {};
for who_idx = 1:metadata.N
    who = metadata.subjects{who_idx};

    for run = 1:metadata.runsPerSubject
        which = which_trials & data.isTrain & data.runId == run & strcmp(data.participant, who);
        
        b = betas(which, :);
        assert(size(b, 1) == metadata.trainingTrialsPerRun);
        
        conds = data.condition(which);
        
        X = [X; mean(b)];
        Y = [Y; conds{1}];
    end
end

t = templateSVM('Standardize',1);
Mdl = fitcecoc(X,Y,'Learners',t,...
    'ClassNames', metadata.contextRoles);

isLoss = resubLoss(Mdl);
%fprintf('          training accuracy: %.3f\n', 1 - isLoss); WRONG -- this is not accuracy!

CVMdl = crossval(Mdl);
oosLoss = kfoldLoss(CVMdl);

%fprintf('          CV accuracy: %.3f\n', 1 - isLoss); WRONG -- this is not accuracy!

% plot confusion matrix
oofLabel = kfoldPredict(CVMdl);
ConfMat = confusionmat(Y,oofLabel);

% Convert the integer label vector to a class-identifier matrix.
[n, p] = size(X);
isLabels = metadata.contextRoles';
nLabels = numel(isLabels)
[~,grpOOF] = ismember(oofLabel,isLabels); 
oofLabelMat = zeros(nLabels,n); 
idxLinear = sub2ind([nLabels n],grpOOF,(1:n)'); 
oofLabelMat(idxLinear) = 1; % Flags the row corresponding to the class 
[~,grpY] = ismember(Y,isLabels); 
YMat = zeros(nLabels,n); 
idxLinearY = sub2ind([nLabels n],grpY,(1:n)'); 
YMat(idxLinearY) = 1; 

figure;
plotconfusion(YMat,oofLabelMat);
h = gca;



%% classifier: multi-class SVM  ... okay per-trial classifier is completely useless
% with leave-run-out cross-validation
% error = 0.5731 w t f
%
accTrain = [];
accTest = [];

for run = 1:metadata.runsPerSubject
    fprintf('------- leave out run = %d\n', run);
    
    which = which_trials & data.isTrain & data.runId ~= run;
    X = betas(which, :);
    Y = data.condition(which);

    Mdl = fitcecoc(X,Y);

    isLoss = resubLoss(Mdl);
%    fprintf('          training accuracy: %.3f\n', 1 - isLoss); WRONG!
%    THIS IS NOT ACCURACY!
    accTrain = [accTrain, 1 - isLoss];
    
    
    % validate on left out run
    %
    which = which_trials & data.isTrain & data.runId == run;
    X = betas(which, :);
    Y = data.condition(which);
    
    y = predict(Mdl, X);
    acc = mean(strcmp(Y, y));
   
%    fprintf('          test accuracy: %.3f\n', acc); WRONG!! THIS IS NOT
%    ACCURACY!
    accTest = [accTest, acc];
end


%% neural network classifier
% accuracy = chance... per-trial tho!
%

inputs = betas;
inputs = score(:,1); % or PC1 instead

targets = nan(size(inputs, 1), 3);
for i = 1:size(inputs, 1)
    switch data.condition{i}
        case 'irrelevant'
            targets(i,:) = [1 0 0];
        case 'modulatory'
            targets(i,:) = [0 1 0];
        case 'additive'
            targets(i,:) = [0 0 1];
    end
end

% only useful trials
%inputs = inputs(which_trials & data.isTrain, :);
%targets = targets(which_trials & data.isTrain, :);


% nevermind -- average the betas
% and not even betas but PC1
%{
X = [];
Y = [];
for who_idx = 1:metadata.N
    who = metadata.subjects{who_idx};

    for run = 1:metadata.runsPerSubject
        which = which_trials & data.isTrain & data.runId == run & strcmp(data.participant, who);
        
        b = betas(which, :);
        assert(size(b, 1) == metadata.trainingTrialsPerRun);
        
        tar = targets(which,:);
        
        %X = [X; mean(b)];
        X = [X; mean(score(which, 1))];
        Y = [Y; tar(1,:)];
    end
end
inputs = X;
targets = Y;
%}


% sanity check
[~, id] = max(targets, [], 2);
anovan(inputs, id);

figure;
hold on;
for i = 1:3
    subplot(1,3,i);
    hist(inputs(id == i), 30);
    title(metadata.contextRoles{i});
end
%{
H = [inputs(id == 1), inputs(id == 2), inputs(id == 3)];
[ns,xs] = hist(H);
hold on;
bar(xs, ns(:,1), 'facecolor', 'blue', 'facealpha', 0.4, 'edgecolor', 'none');
bar(xs, ns(:,2), 'facecolor', 'yellow', 'facealpha', 0.4, 'edgecolor', 'none');
hold off;
box off;
%}

% retarded classifier
guess = nan(size(id));
guess(inputs > 0) = 2;
guess(inputs < 0) = 1;
fprintf('dumb classifier accuracy: %f\n', sum(guess == id) / numel(guess));





inputs = inputs'; % ugh MATLAB
targets = targets';

% from https://github.com/tomov/food-recognition/blob/master/neural_train.m

% Create a Pattern Recognition Network
hiddenLayerSize = 8; % TODO param
net = patternnet(hiddenLayerSize);

% Set up Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
net.trainParam.showWindow = false; % don't show GUI on NCF

% Train the Network
[net,tr] = train(net,inputs,targets);

% Test the Network
outputs = net(inputs);
errors = gsubtract(targets,outputs);
performance = perform(net,targets,outputs);

% View the Network
%view(net)

% View confusion matrix
[c,cm,ind,per] = confusion(targets,outputs);

% patternnet wants column feature vectors, so we rotated those
% but the rest of the code expects them to be rotated
% so we undo that...
% 
inputs = inputs';
targets = targets';
outputs = outputs';

accuracy = classify_get_accuracy(outputs, targets);
fprintf('Success rate = %.2f%%\n', accuracy);    

classifier = net;





%% classify w/ multi-class GLM
% jk doesn't work b/c of macOS Sierra... but i thought I fixed this?? well
% sheeeeiiit...
%

inputs = betas;

targets = nan(size(inputs, 1), 3);
for i = 1:size(inputs, 1)
    switch data.condition{i}
        case 'irrelevant'
            targets(i,:) = [1 0 0];
        case 'modulatory'
            targets(i,:) = [0 1 0];
        case 'additive'
            targets(i,:) = [0 0 1];
    end
end

% use important trials only TODO crossvalidate
%training_trials = which_trials & data.isTrain & data.runId ~= 9;
%test_trails = which_trials & data.isTrain & data.runId == 9;

inputs = inputs(which_trials & data.isTrain, :);
targets = targets(which_trials & data.isTrain, :);

% from classfy_train.m
opts.alpha = 1; % 0 = ridge penalty; 1 = lasso penalty (force betas to zero); default = 1
opts.mtype = 'ungrouped'; % 'grouped' = for multinomial, all betas are in our out together; default = 'ungrouped'
opts.nlambda = 1000; % # lambdas to use; default = 100
opts.lambda_min = 0.00000001; % as a fraction of lambda_max which is derived from the data; default = 0.0001
options = glmnetSet(opts);

% each run is a separate fold TODO rm
%
foldid = data.runId(which_trials & data.isTrain);
%foldid = arrayfun(@(x) find(x == runs), foldid);
disp('folds:');
disp(foldid);

% x = inputs
% y = targets
%
parallel = false;
keep = true;
CVfit = cvglmnet(inputs, targets, 'multinomial', options, 'deviance', [], foldid, parallel, keep);
disp(CVfit);

outputs = cvglmnetPredict(CVfit, inputs, CVfit.lambda_1se, 'response');

accuracy = classify_get_accuracy(outputs, targets);
fprintf('Success rate (lambda = %.4f) is %.2f%%\n', CVfit.lambda_1se, accuracy);

classifier = CVfit;


%{






%% alternative -- run PCA for each subject separately
% CURIOUSLY -- almost the same result!
%
score = [];

figure;

who_idx = 0;
for who = metadata.subjects
    who_idx = who_idx + 1;
    which = which_trials & strcmp(data.participant, who);
    [coeff,s,latent,tsquared,explained,mu] = pca(betas(which, :)); 
    
    score(which, :) = s;
    
    subplot(1, metadata.N, who_idx);
    plot(explained, '-o');
    if who_idx == 1
        ylabel('% variance explained');
    end
end


%% plot PC1 over time for each condition, for each subject for each block
%

for condition = metadata.contextRoles

    figure;

    title(condition{1});
    hold on;
    
    who_idx = 0;
    for who = metadata.subjects
        who_idx = who_idx + 1;

        run_idx = 0;
        for run = 1:metadata.runsPerSubject
            which = which_trials & data.isTrain & data.runId == run & strcmp(data.contextRole, condition) & strcmp(data.participant, who);
            if sum(which) == 0, continue; end % only blocks in that condition
            run_idx = run_idx + 1;
            %assert(sum(which) == metadata.trainingTrialsPerRun);

            subplot(metadata.runsPerContext, metadata.N, (run_idx - 1) * metadata.N + who_idx);

            plot(score(which, 1), '.-');

            set(gca, 'xtick', []);
            set(gca, 'ytick', []);
            if run_idx == 1
                title(who{1});
            end
        end

    end

    hold off;
end


%% plot top two PCs for each block for each subject, separted by condition, iteratively (one trial at a time -- press space)
%

for condition = metadata.contextRoles

    figure;

    title(condition{1});
    hold on;
    
    for t = 1:metadata.trainingTrialsPerRun

        who_idx = 0;
        for who = metadata.subjects
            who_idx = who_idx + 1;

            run_idx = 0;
            for run = 1:metadata.runsPerSubject
                which = which_trials & data.isTrain & data.trialId <= t & data.runId == run & strcmp(data.contextRole, condition) & strcmp(data.participant, who);
                if sum(which) == 0, continue; end % only blocks in that condition
                run_idx = run_idx + 1;
                %assert(sum(which) == metadata.trainingTrialsPerRun);

                subplot(metadata.runsPerContext, metadata.N, (run_idx - 1) * metadata.N + who_idx);

                plot(score(which, 1), score(which, 2), '.-');
                
                set(gca, 'xtick', []);
                set(gca, 'ytick', []);
                if run_idx == 1
                    title(who{1});
                end
%                break;
            end

            %break;
        end
        
        pause;
    end

    hold off;
end


%}
