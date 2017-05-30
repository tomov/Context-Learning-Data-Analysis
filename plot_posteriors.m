function plot_posteriors(data, metadata, params, which_structures)

% Plot the posteriors of the model given subject data.
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

if nargin < 3 || isempty(params)
    params = [0.1249 2.0064]; % by default, use the params from the pilot fit
end

% First simulate the subjects with the causal structure model
%
simulated = simulate_subjects(data, metadata, params, which_structures);

next_subplot_idx = 1; % so you can reorder them by simply rearranging the code

s_id = 0;
for who = metadata.subjects
    s_id = s_id + 1;

    for run = 1:metadata.runsPerSubject

        which = data.which_rows & data.isTrain & strcmp(who, data.participant) & data.runId == run;
        P = simulated.P(data.which_rows & data.isTrain & strcmp(who, data.participant) & data.runId == run, logical(which_structures));
        
        condition = data.contextRole(which);
        condition = condition{1};

        subplot(metadata.N, metadata.runsPerSubject, next_subplot_idx);
        next_subplot_idx = next_subplot_idx + 1;

        plot(P, '-', 'LineWidth', 2);
        text(6, 0.5, condition);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);

        Ms = {}; % only include models we care about
        if which_structures(1), Ms = [Ms, {'M1'}]; end
        if which_structures(2), Ms = [Ms, {'M2'}]; end
        if which_structures(3), Ms = [Ms, {'M3'}]; end
        if which_structures(4), Ms = [Ms, {'M4'}]; end
        
        if run == 1
            ylabel(num2str(s_id));
        end
        if strcmp(who{1}, metadata.subjects{end})
            xlabel('trial');
            if run == metadata.runsPerSubject
                legend(Ms);
            end
        end
        if strcmp(who{1}, metadata.subjects{1})
            title(['Run #', num2str(run)]);
        end
    end
end

