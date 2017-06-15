function T = plot_classifierVsPosterior(data, metadata, params, which_structures)

% Compare classifier P(condition | neural data on trial n) and posterior from model
% P(causal strucutre | h_1:n).
% First train a classifier to predict condition based on activity at trial onset
% separately for each subject.
%
%

utils; % include some libs

masks = {fullfile('masks', 'hippocampus.nii'), ...
         fullfile('masks', 'ofc.nii'), ...
         fullfile('masks', 'striatum.nii'), ...
         fullfile('masks', 'vmpfc.nii'), ...
         fullfile('masks', 'bg.nii'), ...
         fullfile('masks', 'pallidum.nii'), ...
         fullfile('masks', 'visual.nii'), ...
         fullfile('masks', 'v1.nii'), ...
         fullfile('masks', 'motor.nii'), ...
         fullfile('masks', 'sensory.nii')};
z_score = 'z-run';

% notice the structures are hardcoded here -- this is b/c we always classify 3 conditions
%
assert(isequal(which_structures, [1 1 1 0]));
which_structures = logical(which_structures);
structure_names = {'M1', 'M2', 'M3', 'M4'};

% Simulate the subjects
%
simulated = simulate_subjects(data, metadata, params, which_structures);

% results table
%
table_rowNames = {};
table_colNames = {'acc', 'MCVE', 'KL_P_Q', 'KL_U_Q', 'KL_p', 'Bhat_P_Q', 'Bhat_U_Q', 'Bhat_p','Hell_P_Q', 'Hell_U_Q', 'Hell_p'};
table = [];

next_subplot_idx = 1;

% we're comparing multiple masks
%
for mask = masks
    mask = mask{1};
    [~, maskname, ~] = fileparts(mask);

    % Load classifiers for all subjects, the targets and outputs. Note that
    % each row in all_targets and all_outputs corresponds to a row in data (i.e. all
    % trials are here, not just a subset of them, so we can refer to rows in all_targets
    % and all_outputs using logical vectors from data)
    %
    [classifiers, all_targets, all_outputs] = classify_each_subject(mask, z_score);

    % We assume that all good subjects are included an no others
    %
    which_subjects = ~cellfun(@isempty, classifiers);
    assert(isequal(metadata.allSubjects(which_subjects), metadata.subjects));
    which_rows = data.which_rows & data.isTrain; % only look at training trials
    subjects = find(which_subjects);
    
    %
    % Plots
    %
    
    % Classifier prediction P(condition | neural data on trial n)
    %
    for condition = metadata.contextRoles
        condition = condition{1};

        P_outputs = [];
        for n = 1:metadata.trainingTrialsPerRun
            which = which_rows & data.isTrain & strcmp(data.contextRole, condition) & data.trialId == n;

            P_outputs_n = all_outputs(which, which_structures);
            P_outputs_n = mean(P_outputs_n, 1);
            P_outputs = [P_outputs; P_outputs_n];
        end

        subplot(numel(masks) + 1, 3, next_subplot_idx);
        next_subplot_idx = next_subplot_idx + 1;
        plot(P_outputs, 'o-', 'LineWidth', 2);
        if isequal(condition, metadata.contextRoles{1})
            ylabel(maskname);
        end
        if isequal(mask, masks{1})
            title([condition, ' runs']);
        end
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        
        if isequal(condition, metadata.contextRoles{end})
            legend(metadata.contextRoles{which_structures});
        end
    end
        
    %
    % Compute statistics
    %
    
    % Compute accuracy
    %
    accuracy = classify_get_accuracy(all_outputs(which_rows, :), all_targets(which_rows, :));
    fprintf('Success rate for %s is %.2f%%\n', maskname, accuracy);

    % Compute mean cross-validated error for lambda within 1 standard error
    % of the minimum lambda, averaged across subjects
    %
    MCVEs = [];
    for subj_id = subjects
        c = classifiers{subj_id};
        MCVE = c.cvm(find(c.lambda == c.lambda_1se));
        MCVEs = [MCVEs, MCVE];
    end
    accuracy = classify_get_accuracy(all_outputs(which_rows, :), all_targets(which_rows, :));
    fprintf('mean mean cross-validated error %s is %.2f\n', maskname, mean(MCVEs));

    
    % Distributions we're comparing
    %
    P = simulated.P(which_rows, which_structures); % posterior from model
    Q = all_outputs(which_rows, :); % condition PMF from classifier
    unif = ones(size(P)) / size(P, 2); % uniform
    
    % Compute KL divergence
    %
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
    
    %TODO try triangular discrimination 
    % http://www.sciencedirect.com/science/article/pii/S1053811910003915
    
    row = [accuracy, mean(MCVEs)...
           mean(KL_P_Q), mean(KL_unif_Q), KL_p, ...
           mean(Bhat_P_Q), mean(Bhat_Q_unif), Bhat_p, ...
           mean(Hell_P_Q), mean(Hell_Q_unif), Hell_p];
    table = [table; row];
    table_rowNames = [table_rowNames, {maskname}];
    
    save(fullfile('classifier', ['classifierVsPosterior_', maskname, '_', z_score, '.mat']));
end

T = array2table(table, 'RowNames', table_rowNames, 'VariableNames', table_colNames);


% Plot model posteriors P(M | h_1:n) for comparison
%
for condition = metadata.contextRoles
    condition = condition{1};

    P = [];
    for n = 1:metadata.trainingTrialsPerRun
        which = which_rows & data.isTrain & strcmp(data.contextRole, condition) & data.trialId == n;

        P_n = simulated.P(which, which_structures);
        P_n = mean(P_n, 1);
        P = [P; P_n];
    end

    subplot(numel(masks) + 1, 3, next_subplot_idx);
    next_subplot_idx = next_subplot_idx + 1;
    plot(P, 'o-', 'LineWidth', 2);
    xlabel('n (trial #)');
    if isequal(condition, metadata.contextRoles{1})
        ylabel('P(M | h_{1:n})');
    end
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);

    if isequal(condition, metadata.contextRoles{end})
        legend(structure_names{which_structures});
    end
end

