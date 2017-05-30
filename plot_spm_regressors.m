function plot_spm_regressors(EXPT, glmodel, subj, run)

% Plot the regressors for a given glmodel, subject & run.
% Requires the SPM structure to have been generated i.e. the single-subject
% GLM should have been run
%

include_motion = false;
%mask = 'hippocampus.nii';

multi = EXPT.create_multi(glmodel, subj, run);
load(fullfile('temp', 'context_create_multi.mat')); % WARNING WARNING WARNING: MASSIVE COUPLING. This relies on context_create_multi saving its state into this file. I just don't wanna copy-paste or abstract away the code that load the data from there
disp(condition)

modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
load(fullfile(modeldir,'SPM.mat'));

TR = EXPT.TR;
sess_prefix = ['Sn(', num2str(run), ')'];
trs = 1 : TR : TR*length(SPM.Sess(run).row); % or start at 0? how does slice timing interpolation work in SPM?


figure;

% which regressors to display
%
cols = SPM.Sess(run).col;
if ~include_motion
    cols = cols(1:end-6); % ditch motion regressors
end

% show stimulus sequence
% this is ConLearn-specific 
%
subplot(length(cols) + 1, 1, 1);
hold on;

is_sick = strcmp(data.sick(which_rows), 'Yes');
is_train = which_train(which_rows);
resps = data.response.keys(which_rows);
chose_sick = double(strcmp(resps, 'left'));
is_timeout = strcmp(resps, 'None');
corr = data.response.corr(which_rows);

onsets = cellfun(@str2num, data.actualChoiceOnset(which_rows)');
for t=onsets
    plot([t t], [-2 2], '--', 'Color', [0.8 0.8 0.8]);
end

feedback_onsets = cellfun(@str2num, data.actualFeedbackOnset(which_train)');
for i=1:length(feedback_onsets)
	if corr(i)
        color = 'blue';
    else
        color = 'red';
    end
    plot([feedback_onsets(i) feedback_onsets(i)], [-2 2], 'Color', color);
end

stims = strcat('x', num2str(data.cueId(which_rows) + 1), 'c', num2str(data.contextId(which_rows) + 1));
show_text_nontimeout = @(which_trials, color) text(onsets(which_trials & ~is_timeout), chose_sick(which_trials & ~is_timeout), stims(which_trials & ~is_timeout,:), 'Color', color);
show_text_timeout = @(which_trials, color) text(onsets(which_trials & is_timeout), 0.5 * ones(sum(which_trials & is_timeout), 1), stims(which_trials & is_timeout,:), 'Color', color);

% sick training trials
show_text_nontimeout(is_sick & is_train, [0 0.5 0]);
show_text_timeout(is_sick & is_train, [0 0.5 0]);
% not sick training trials
show_text_nontimeout(~is_sick & is_train, [0.5 0.5 0]);
show_text_timeout(~is_sick & is_train, [0.5 0.5 0]);
% test trials
show_text_nontimeout(~is_train, [0.3 0.3 0.3]);
show_text_timeout(~is_train, [0.3 0.3 0.3]);

ylim([-0.3 1.1]);
yticks([0 0.5 1]);
yticklabels({'chose not sick', 'timeout', 'chose sick'});

% for legend
h = [
     plot(NaN,NaN, '--', 'Color', [0.8 0.8 0.8],'LineWidth', 1); ...
     plot(NaN,NaN, 'Color', 'red','LineWidth', 1); ...
     plot(NaN,NaN, 'Color', 'blue','LineWidth', 1); ...
     plot(NaN,NaN, 'Color', [0 0.5 0],'LineWidth', 2); ...
     plot(NaN,NaN,'Color', [0.5 0.5 0],'LineWidth', 2); ...
     plot(NaN,NaN,'Color', [0.3 0.3 0.3],'LineWidth', 2)];
legend(h, {'trial onset', 'WRONG/TIMEOUT feedback', 'CORRECT feedback', 'sick outcome', 'not sick outcome', 'test trial'});

title(sprintf('Model %d, subject %d, run %d, condition: %s', glmodel, subj, run, upper(condition)));
hold off;


% iterate over regressors
%
plot_idx = 1;
for i = cols
    assert(strncmp(SPM.xX.name{i}, sess_prefix, length(sess_prefix)));
    
    plot_idx = plot_idx + 1;
    subplot(length(cols) + 1, 1, plot_idx);
    hold on;
    h = [];
    
    % plot trial onsets / offsets as vertical dashed lines
    %
    feedback_onsets = cellfun(@str2num, data.actualChoiceOnset(which_rows)');
    for t=feedback_onsets
        plot([t t], [-1 1], '--', 'Color', [0.8 0.8 0.8]);
    end
    
    % plot original regressor from model
    %
    eps = 1e-6;
    leg = {}; % legend
    for j = 1:length(multi.names)
        n = length(multi.onsets{j});
        onsets = multi.onsets{j};
        durations = multi.durations{j};
        if ~(size(durations, 1) == size(onsets, 1) && size(durations, 2) == size(onsets, 2))
            % need to rotate one of the vectors
            durations = durations';
        end
        assert(size(durations, 1) == size(onsets, 1) && size(durations, 2) == size(onsets, 2));
        x = [onsets - eps; onsets; onsets + durations; onsets + durations + eps];
        x = x(:)';
        x = [0 x max(trs)];
        if ~isempty(strfind(SPM.xX.name{i}, [' ', multi.names{j}, '*']))
            y = [zeros(1,n); ones(1,n); ones(1,n); zeros(1,n)];
            y = y(:)';
            y = [0 y 0];
            h = [h, plot(x, y, 'LineWidth', 1, 'Color', 'red')];
            leg = [leg; multi.names{j}];
        end
            
        if isfield(multi, 'pmod') && j <= length(multi.pmod)
            for k = 1:length(multi.pmod(j).name)
                if ~isempty(strfind(SPM.xX.name{i}, ['x', multi.pmod(j).name{k}, '^']))
                    y = reshape(multi.pmod(j).param{k}, 1, n);
                    y = [zeros(1,n); y; y; zeros(1,n)];
                    y = y(:)';
                    y = [0 y 0];
                    h = [h, plot(x, y, 'LineWidth', 1, 'Color', 'red')];
                    leg = [leg; ['pmod: ', multi.pmod(j).name{k}]];
                end
            end
        end
    end
    
    % plot regressor convolved with HRF
    %
    h = [h, plot(trs, SPM.xX.X(SPM.Sess(run).row, i)', 'Color', 'blue')];
    leg = [leg; {SPM.xX.name{i}}];
    
    % TODO plot beta
    %
    %beta_vec = ccnl_get_beta(EXPT, glmodel, i, mask, [subj]);
    
    yL = get(gca,'YLim');
    ylim([yL(1), yL(2) + 0.1]);    
    title(SPM.xX.name{i}, 'Interpreter', 'none');
    legend(h, leg, 'Interpreter', 'none');
        
    hold off
end
