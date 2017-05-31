function [classifiers, all_targets, all_outputs] = classify_single_subject(mask, z_score)

% Train a separate classifier for each subject to predict the run condition
% based on neural activation.
%

method = 'cvglmnet';
%z_score = 'z-run'; -- passed as argument
predict_what = 'condition';
runs = [1:9];
trials = [1:24];

% assumes we're working with all good subjects
%
subjects = getGoodSubjects();

% Run computation for each mask
%
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
 
