function plot_classifierPredictions(data, metadata, mask)

% Plot classifier P(condition | neural data on trial n) for each run for
% each subject.
%

z_score = 'z-run';

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

next_subplot_idx = 1; % so you can reorder them by simply rearranging the code

% Get classifier predictions
%
s_ord = 0;
for subj = subjects
    s_ord = s_ord + 1;
    who = metadata.allSubjects{subj};

    for run = 1:metadata.runsPerSubject

        which = which_rows & data.isTrain & strcmp(who, data.participant) & data.runId == run;        

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
                legend(metadata.contextRoles{logical([1 1 1 0])});
            end
        end
        if strcmp(who, metadata.subjects{1})
            title(['Run #', num2str(run)]);
        end
    end
end
