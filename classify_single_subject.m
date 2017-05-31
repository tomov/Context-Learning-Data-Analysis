%function plot_classifierVsPosterior(data, metadata, params, which_structures)

% Compare classifier P(condition | neural data on trial n) and posterior from model
% P(causal strucutre | h_1:n).
% First train a classifier to predict condition based on activity at trial onset
% separately for each subject.
% Based on classify_single_subject.sh which overwhelmed NCF lol
%
%

utils; % include some libs

% notice the structures are hardcoded here -- this is b/c we always classify 3 conditions
%
assert(isequal(which_structures, [1 1 1 0]));
which_structures = logical(which_structures);
structure_names = {'M1', 'M2', 'M3', 'M4'};

% Simulate the subjects
%
simulated = simulate_subjects(data, metadata, params, which_structures);


method = 'cvglmnet';
masks = {fullfile('masks', 'hippocampus.nii'), ...
         fullfile('masks', 'ofc.nii'), ...
         fullfile('masks', 'striatum.nii'), ...
         fullfile('masks', 'vmpfc.nii'), ...
         fullfile('masks', 'bg.nii'), ...
         fullfile('masks', 'pallidum.nii'), ...
         fullfile('masks', 'visual.nii'), ...
         fullfile('masks', 'motor.nii'), ...
         fullfile('masks', 'sensory.nii')};
%masks = {masks{1}};
z_score = 'z-run';
predict_what = 'condition';
runs = [1:9];
trials = [1:24];

% assumes we're working with all good subjects
%
subjects = getGoodSubjects();
assert(isequal(metadata.allSubjects(subjects), metadata.subjects));

% results table
%
table_rowNames = {};
table_colNames = {'acc', 'KL_P_Q', 'KL_U_Q', 'KL_p', 'Bhat_P_Q', 'Bhat_U_Q', 'Bhat_p','Hell_P_Q', 'Hell_U_Q', 'Hell_p'};
table = [];

% Run computation for each mask
%
for mask = masks
    mask = mask{1};
    [~, maskname, ~] = fileparts(mask);
    
    fprintf('\n\n----- Mask %s ---\n\n', maskname);
    
    % Load the classifiers, if they have been trained already
    % Or train them anew and save them if they haven't
    %
    classifiers_filename = fullfile('classifier', ['single_subject_', method, '_', z_score, '_', predict_what, '_', maskname, '.mat']);
    fprintf('filename = %s\n', classifiers_filename);
    
    if exist(classifiers_filename, 'file') ~= 2
        
        classifiers = {}; % classifier for each subject
        all_targets = []; % classifier targets for all subjects / rows of the data
        all_outputs = []; % classifier outputs for all subjects / rows of the data
        for subj = subjects
            fprintf('  ---- Classifying subject %d ---------\n', subj);

            [classifier, ~, targets, outputs, which_rows] = classify_train(method, runs, trials, subj, mask, predict_what, z_score);
            all_targets(which_rows, :) = targets;
            all_outputs(which_rows, :) = outputs;
            classifiers{subj} = classifier;
        end

        fprintf('Saving classifiers to %s\n', classifiers_filename);
        save(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets', '-v7.3');
    else
        fprintf('Loading saved classifiers\n');
        load(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets');
    end


    all_targets = []; % classifier targets for all subjects / rows of the data
    all_outputs = []; % classifier outputs for all subjects / rows of the data
    for subj = subjects
        fprintf('  ---- Testing subject %d ---------\n', subj);

        [~, targets, outputs, which_rows] = classify_test(method, classifiers{subj}, runs, trials, subj, mask, predict_what, z_score);
        all_targets(which_rows, :) = targets;
        all_outputs(which_rows, :) = outputs;
    end
    save(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets', '-v7.3');
    
    
    which_rows = data.which_rows & data.isTrain; % only look at training trials
    
    %
    % Plot classifier P(condition | neural data on trial n) vs. P(M | h_1:n) from
    % model
    %
    
    %{
    
    figure;
    next_subplot_idx = 1;
    
    % Classifier predictions
    %
    for condition = metadata.contextRoles

        P_outputs = [];
        for n = 1:metadata.trainingTrialsPerRun
            which = which_rows & data.isTrain & strcmp(data.contextRole, condition) & data.trialId == n;

            P_outputs_n = all_outputs(which, which_structures);
            P_outputs_n = mean(P_outputs_n, 1);
            P_outputs = [P_outputs; P_outputs_n];
        end

        subplot(2, 3, next_subplot_idx);
        next_subplot_idx = next_subplot_idx + 1;
        plot(P_outputs, 'o-', 'LineWidth', 2);
        xlabel('n (trial #)');
        ylabel('P(condition | neural data)');
        title(strcat('Classifier predictions in ', {' '}, condition));

        legend(metadata.contextRoles{which_structures});
    end
    
    % Model posteriors
    %
    for condition = metadata.contextRoles

        P = [];
        for n = 1:metadata.trainingTrialsPerRun
            which = which_rows & data.isTrain & strcmp(data.contextRole, condition) & data.trialId == n;

            P_n = simulated.P(which, which_structures);
            P_n = mean(P_n, 1);
            P = [P; P_n];
        end

        subplot(2, 3, next_subplot_idx);
        next_subplot_idx = next_subplot_idx + 1;
        plot(P, 'o-', 'LineWidth', 2);
        xlabel('n (trial #)');
        ylabel('P(M | h_{1:n})');
        title(strcat('Posterior after each trial for ', {' '}, condition));

        legend(structure_names{which_structures});
    end
    
    title(maskname);
    %}
    
    
    % Compute accuracy
    %
    accuracy = classify_get_accuracy(all_outputs(which_rows, :), all_targets(which_rows, :));
    fprintf('Success rate for %s (lambda = %.4f) is %.2f%%\n', maskname, classifier.lambda_1se, accuracy);

    % Compute KL divergence
    %
    P = simulated.P(which_rows, which_structures); % posterior from model
    Q = all_outputs(which_rows, :); % condition PMF from classifier
    unif = ones(size(P)) / size(P, 2); % uniform
    
    % mean absolute deviation ? or divergence? TODO
    
    KL_P_Q = KL_divergence(P, Q);
    KL_P_unif = KL_divergence(P, unif);
    KL_unif_Q = KL_divergence(unif, Q);
    [~, KL_p] = ttest2(KL_P_Q, KL_unif_Q);
    
    % these are -Infs
    KL_Q_P = KL_divergence(Q, P);
    KL_Q_unif = KL_divergence(Q, unif);
    KL_unif_P = KL_divergence(unif, P);
    
    Bhat_P_Q = Bhattacharyya_distance(P, Q);
    Bhat_P_unif = Bhattacharyya_distance(P, unif);
    Bhat_Q_unif = Bhattacharyya_distance(Q, unif);
    [~, Bhat_p] = ttest2(Bhat_P_Q, Bhat_Q_unif);
    
    Hell_P_Q = Hellinger_distance(P, Q);
    Hell_P_unif = Hellinger_distance(P, unif);
    Hell_Q_unif = Hellinger_distance(Q, unif);
    [~, Hell_p] = ttest2(Hell_P_Q, Hell_Q_unif);
    
    row = [accuracy, ...
           mean(KL_P_Q), mean(KL_unif_Q), KL_p, ...
           mean(Bhat_P_Q), mean(Bhat_Q_unif), Bhat_p, ...
           mean(Hell_P_Q), mean(Hell_Q_unif), Hell_p];
    table = [table; row];
    table_rowNames = [table_rowNames, {maskname}];
    
    
    %{
    
    
    
    figure;
    next_subplot_idx = 1; % so you can reorder them by simply rearranging the code
    
    % Get classifier predictions
    %
    s_ord = 0;
    for subj = subjects
        s_ord = s_ord + 1;
        who = metadata.allSubjects{subj};

        classifier = classifiers{subj};
        
        %accuracy = classify_get_accuracy(outputs, targets);
        %[r, p] = corrcoef(outputs, simulated.P(which_rows, which_structures));
        %fprintf('Success rate (lambda = %.4f) is %.2f%%\n', classifier.lambda_1se, accuracy);
        
        for run = 1:metadata.runsPerSubject

            which = data.which_rows & data.isTrain & strcmp(who, data.participant) & data.runId == run;        
        
            condition = data.contextRole(which);
            condition = condition{1};

            subplot(metadata.N, metadata.runsPerSubject, next_subplot_idx);
            next_subplot_idx = next_subplot_idx + 1;

            plot(all_outputs(which, :), '-', 'LineWidth', 2);
            text(6, 0.5, condition);
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);

            if run == 1
                ylabel(num2str(s_ord));
            end
            if strcmp(who, metadata.subjects{end})
                xlabel('trial');
                if run == metadata.runsPerSubject
                    legend(structure_names{which_structures});
                end
            end
            if strcmp(who, metadata.subjects{1})
                title(['Run #', num2str(run)]);
            end
        end
    end
    
    %}
end
 
T = array2table(table, 'RowNames', table_rowNames, 'VariableNames', table_colNames);