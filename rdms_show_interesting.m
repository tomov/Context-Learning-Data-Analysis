function rdms_show_interesting(rows, cols, all_subject_rhos, all_vs_all, control_model_idxs)


%% Interesting visualizations
% show select neural / model pairs that I think are interesting, like bar
% plots
%
assert(~all_vs_all); % only works for neural vs. model comparisons

neural_idxs = [1, 5, 2, 3, 4]; % ROIs to plot at trial onset
neural_idxs = [neural_idxs; neural_idxs + numel(rows)/2]; % at feedback too
neural_idxs = neural_idxs(:)';

model_idxs = [1:7, 18, 20, 26, 29, 30, 31, 32, 15, 35, 8];

relative_to_time = false; % #KNOB whether to do the analysis/plots relative to time
control_for_time = ~isempty(find(control_model_idxs == 8));
assert(~relative_to_time || ~control_for_time);

all_subject_zs = atanh(all_subject_rhos); % Fisher z transform the rhos

figure;

for i = 1:numel(neural_idxs)
    neural_idx = neural_idxs(i);
    
    % Get fisher z-transformed correlation coefficients
    %
    zs = squeeze(all_subject_zs(neural_idx, model_idxs, :));
    
    if relative_to_time
        time_zs = squeeze(all_subject_zs(neural_idx, 8, :));
        
        % One-sample t-test against the time zs
        %
        [h, ps, ci, stats] = ttest(zs' - time_zs);
        assert(numel(h) == numel(model_idxs));
        
        % plot the z's relative to time
        %
        zs = time_zs - zs;
    else
        % One-sample t-test them against 0
        %
        [h, ps, ci, stats] = ttest(zs');
        assert(numel(h) == numel(model_idxs));
    end
    
    % Compute SEs
    %
    %[sems, means] = wse(zs'); % WSE's
    means = mean(zs, 2);
    sems = (std(zs, [], 2) / sqrt(size(zs, 2)))'; % ...or just normal SEMs
    assert(numel(sems) == numel(model_idxs));
    
    % Plot
    %
    subplot(numel(neural_idxs) + 1, 1, i);
    
    h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    xs = h(1).XData;
    hold on;
    errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    hold off;
    
    set(gca, 'xtick', xs);
    models = {cols(model_idxs).name};
    labels = {};
    for j = 1:numel(model_idxs)
        labels{j} = sprintf('p < 10^{%d}', round(log10(ps(j))));
    end
    xticklabels(labels);
 
    if ~control_for_time
        set(gca, 'ylim', [0 0.405]);
    else
        set(gca, 'ylim', [-0.1 0.03]);
    end
    set(gca, 'ytick', mean(ylim));
    yticklabels({rows(neural_idx).name});
    %ytickangle(30);
    if i == 1
        title('Fisher z-transformed Spearman rank correlation');
    end
    if i == numel(neural_idxs)
        subplot(numel(neural_idxs) + 1, 1, i + 1);
        h = bar(means * 0, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
        set(gca, 'xtick', xs);
        xticklabels(models);
        xtickangle(45);
    end
end


