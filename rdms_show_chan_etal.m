function rdms_show_chan_etal(rows, cols, all_subject_rhos, all_vs_all, control_model_idxs)

%% Chan et al.-like bar plots for their models
%
assert(~all_vs_all); % only works for neural vs. model comparisons

neural_idxs = [1:(numel(rows) / 2)]; % ROIs to plot -- trial_onset
neural_idxs = [1:(numel(rows) / 2)] + (numel(rows) / 2); % ROIs to plot -- feedbac_onset
model_idxs = [1:7, 28, 9, 10, 11, 8]; % models to plot -- Chan et al.
relative_to_time = false; % #KNOB whether to do the analysis/plots relative to time
control_for_time = ~isempty(find(control_model_idxs == 8));
assert(~relative_to_time || ~control_for_time);

all_subject_zs = atanh(all_subject_rhos); % Fisher z transform the rhos

figure;

for i = 1:numel(model_idxs)
    model_idx = model_idxs(i);
    
    % Get fisher z-transformed correlation coefficients
    %
    zs = squeeze(all_subject_zs(neural_idxs, model_idx, :));
    
    if relative_to_time
        time_zs = squeeze(all_subject_zs(neural_idxs, 8, :));
        
        % Paired-sample t-test them against the time zs
        %
        [h, ps, ci, stats] = ttest(zs', time_zs');
        assert(numel(h) == numel(neural_idxs));
        
        % plot the z's relative to time
        %
        zs = time_zs - zs;
    else
        % One-sample t-test them against 0
        %
        [h, ps, ci, stats] = ttest(zs');
        assert(numel(h) == numel(neural_idxs));
    end
    
    % Compute SEs
    %
    %[sems, means] = wse(zs'); % WSE's
    means = mean(zs, 2);
    sems = (std(zs, [], 2) / sqrt(size(zs, 2)))'; % ...or just normal SEMs
    assert(numel(sems) == numel(neural_idxs));
    
    % Plot
    %
    subplot(numel(model_idxs), 1, i);
    
    h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    xs = h(1).XData;
    hold on;
    errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    hold off;
    
    set(gca, 'xtick', xs);
    rois = {rows(neural_idxs).name};
    if i == numel(model_idxs)
        xticklabels(rois);
        xtickangle(45);
    else
        labels = {};
        for j = 1:numel(neural_idxs)
            labels{j} = sprintf('p < 10^{%d}', round(log10(ps(j))));
        end
        xticklabels(labels);
    end
 
    if ~control_for_time
        set(gca, 'ylim', [0 0.405]);
    else
        set(gca, 'ylim', [-0.2 0.05]);
    end
    set(gca, 'ytick', mean(ylim));
    yticklabels({cols(model_idx).name});
    %ytickangle(30);
    if i == 1
        title('Fisher z-transformed Spearman rank correlation');
    end
end

