function plot_PE(data, metadata, params, which_structures)

% Plot the PE of the model given subject data. For all subjects and
% all runs separately.
% 
% INPUT 
% data, metadata = subject data and metadata as output by load_data
% params = vector of hyperparameters:
%    params(:, 1) = prior_variance
%    params(:, 2) = inv_softmax_temp
%          for fixed effects, only have 1 row and these parameters will be
%          used for all subjects. For random effects, have a separate row
%          with the parameters for each subject
% which_structures = which causal structures to use, as a logical vector,
%                    e.g. [1 1 1 0] = M1, M2, M3
%

if nargin < 4 || isempty(which_structures)
    which_structures = [1 1 1 0]; % by defulat, use M1 M2 M3
end

% First simulate the subjects with the causal structure model
%
simulated = simulate_subjects(data, metadata, params, which_structures);

next_subplot_idx = 1; % so you can reorder them by simply rearranging the code

s_ord = 0;
for who = metadata.subjects
    s_ord = s_ord + 1;

    for run = 1:metadata.runsPerSubject

        which = data.which_rows & data.isTrain & strcmp(who, data.participant) & data.runId == run;
        values = simulated.values(which);
        new_values = simulated.new_values(which);
        outcomes = data.outcome(which);
        
        PE = (outcomes' - values').^2; % alternative 1
        % PE = (new_values' - values').^2; % alternative 2
        
        condition = data.contextRole(which);
        condition = condition{1};

        subplot(metadata.N, metadata.runsPerSubject, next_subplot_idx);
        next_subplot_idx = next_subplot_idx + 1;

        plot(PE, '-', 'LineWidth', 2);
        text(6, mean(ylim), condition);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);

        if run == 1
            ylabel(num2str(s_ord));
        end
        if strcmp(who{1}, metadata.subjects{end})
            xlabel('trial');
            if run == metadata.runsPerSubject
                legend({'PE^2'});
            end
        end
        if strcmp(who{1}, metadata.subjects{1})
            title(['Run #', num2str(run)]);
        end
    end
end

