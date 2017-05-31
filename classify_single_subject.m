% Compare classifier P(condition | neural data on trial n) and posterior from model
% P(causal strucutre | h_1:n).
% First train a classifier to predict condition based on activity at trial onset
% separately for each subject.
% Based on classify_single_subject.sh which overwhelmed NCF lol
%

% TODO parametrize and take in data, metadata, params ,etc as inputs 
% just like other shit in plot_figure

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts_100_prior.mat'), 'results', 'results_options');
params = results(1).x;
options = results_options(1);
disp('Using parameters:');
disp(params);
disp('generated with options:');
disp(options);
%assert(options.isFmriData == true);
%assert(~options.fixedEffects);
assert(isequal(options.which_structures, [1 1 1 0]));

% notice the structures are hardcoded here -- this is b/c we always classify 3 conditions
%
which_structures = logical([1 1 1 0]);
structure_names = {'M1', 'M2', 'M3', 'M4'};

% Simulate the subjects
%
simulated = simulate_subjects(data, metadata, params, which_structures);


method = 'cvglmnet';
masks = { ... %fullfile('masks', 'hippocampus.nii'), ...
         fullfile('masks', 'ofc.nii'), ...
         fullfile('masks', 'striatum.nii'), ...
         fullfile('masks', 'vmpfc.nii'), ...
         fullfile('masks', 'bg.nii'), ...
         fullfile('masks', 'pallidum.nii'), ...
         fullfile('masks', 'visual.nii'), ...
         fullfile('masks', 'motor.nii'), ...
         fullfile('masks', 'sensory.nii')};
%masks = {masks{2}};
z_score = 'z-none';
predict_what = 'condition';
runs = [1:9];
trials = [1:24];

subjects = getGoodSubjects();
assert(isequal(metadata.allSubjects(subjects), metadata.subjects)); % for the visualization

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
    
    figure;
    
    next_subplot_idx = 1; % so you can reorder them by simply rearranging the code
    
    % Get classifier predictions
    %
    all_targets = []; % classifier targets for all subjects / rows of the data
    all_outputs = []; % classifier outputs for all subjects / rows of the data
    s_ord = 0;
    for subj = subjects
        s_ord = s_ord + 1;
        who = metadata.allSubjects{subj};

        classifier = classifiers{subj};
        [inputs, targets, outputs, which_rows] = classify_test(method, classifier, runs, trials, subj, mask, predict_what, z_score);
        all_targets(which_rows, :) = targets;
        all_outputs(which_rows, :) = outputs;
        
        %accuracy = classify_get_accuracy(outputs, targets);
        %[r, p] = corrcoef(outputs, simulated.P(which_rows, which_structures));
        %fprintf('Success rate (lambda = %.4f) is %.2f%%\n', classifier.lambda_1se, accuracy);
        
        for run = 1:metadata.runsPerSubject

            which = which_rows & data.isTrain & strcmp(who, data.participant) & data.runId == run;        
            which_outputs = which(which_rows); % outputs only represent rows for that subject
        
            condition = data.contextRole(which_rows);
            condition = condition{1};

            subplot(metadata.N, metadata.runsPerSubject, next_subplot_idx);
            next_subplot_idx = next_subplot_idx + 1;

            plot(outputs(which_outputs, :), '-', 'LineWidth', 2);
            text(6, 0.5, condition);
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);

            if run == 1
                ylabel(num2str(s_ord));
            end
            if strcmp(who, metadata.subjects{end})
                xlabel('trial');
                if run == metadata.runsPerSubject
                    legend(struct_names{which_structures});
                end
            end
            if strcmp(who, metadata.subjects{1})
                title(['Run #', num2str(run)]);
            end
        end
    end
    
    save(classifiers_filename, 'classifiers', 'all_outputs', 'all_targets', '-v7.3');

end
 