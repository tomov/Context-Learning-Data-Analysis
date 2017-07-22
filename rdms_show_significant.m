function rdms_show_significant(rows, cols, table_Rho, table_P, all_subject_rhos, lmes, do_LME, all_vs_all, bonferroni)

% Show the significant positive correlations from the second-order RDM
%
% INPUT:
% rows, cols = struct arrays of RDMs as input to rdms_second_order()
% table_Rho, table_P, all_subject_rhos, do_LME = as output by rdms_second_order()
% all_vs_all = whether we're comparing all RDMs vs. themselves, or Neural vs. Model RDMs
% bonferroni = whether to use bonferroni correction for the significance test
%

if bonferroni
    which = table_Rho > 0 & table_P < 0.05 / numel(table_P); % Bonferroni correction
else
    which = table_Rho > 0 & table_P < 0.01;
end

assert(~all_vs_all); % only works for neural vs. model comparisons

all_subject_zs = atanh(all_subject_rhos); % Fisher z transform the rhos

figure;

plot_idx = 0;
for row_idx = 1:numel(rows)
    % Some stats -- which correlations are significant?
    %
    if do_LME 
        % LMEs
        %
        col_idxs = [];
        means = [];
        sems = [];
        ps = [];
        for col_idx = 1:numel(cols)
            if size(lmes, 1) < row_idx || size(lmes, 2) < col_idx || isempty(lmes{row_idx, col_idx}), continue; end
            lme = lmes{row_idx, col_idx};
            p = lme.Coefficients(2, 'pValue').pValue;
            coef = lme.Coefficients(2, 'Estimate').Estimate;
            if coef <= 0 || p >= 0.05, continue; end
            
            means = [means, coef];
            sems = [sems, (lme.Coefficients(2, 'Upper').Upper - lme.Coefficients(2, 'Lower').Lower) / 2];
            ps = [ps, p];
            col_idxs = [col_idxs, col_idx];
        end
        
        if isempty(col_idxs), continue; end
        
    else
        % Run the t-tests manually with the Fisher z-transformed
        % correlation coefficients
        %
        col_idxs = find(which(row_idx,:))';

        if isempty(col_idxs), continue; end

        % Get fisher z-transformed correlation coefficients
        %
        zs = reshape(all_subject_zs(row_idx, col_idxs, :), [numel(col_idxs) size(all_subject_zs, 3)]);

        % One-sample t-test them against 0
        %
        [h, ps, ci, stats] = ttest(zs');
        assert(numel(h) == numel(col_idxs));

        % Compute SEs
        %
        %[sems, means] = wse(zs'); % WSE's
        means = mean(zs, 2);
        sems = (std(zs, [], 2) / sqrt(size(zs, 2)))'; % normal SEMs
        assert(numel(sems) == numel(col_idxs));
        assert(numel(means) == numel(col_idxs));
    end

    % Plot
    %
    plot_idx = plot_idx + 1;
    if do_LME 
        subplot(4, 7, plot_idx);
    else
        if bonferroni
            subplot(2, 5, plot_idx);
        else
            subplot(3, 7, plot_idx);
        end
    end

    % Plot the bar graphs with error bars
    %
    h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    xs = h(1).XData;
    hold on;
    errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    hold off;
    
    % Put the p-values
    %
    for i = 1:numel(col_idxs)
        significance = @(p) repmat('*', 1, floor(-log10(p)) - 1);
        if bonferroni
            stars = significance(ps(i) * numel(table_P)); % Bonferroni correction
        else
            stars = significance(ps(i));
        end
        %stars = int2str(length(stars));
        text(xs(i) - length(stars) * 0.3 / 2, means(i) + sems(i) * 1.2, stars, 'Rotation', 45);
    end

    % Put the ROI / model names and figure title
    %
    set(gca, 'xtick', xs);
    models = {cols(col_idxs).name};
    labels = {};
    for j = 1:numel(col_idxs)
        labels{j} = sprintf('p < 10^{%d}', round(log10(ps(j))));
    end
    xticklabels(models);
    xtickangle(55);

    t = rows(row_idx).name;
    if t(end) == 'f'
        t = [t(1:end-2), '_{feedback\_onset}'];
    else
        t = [t(1:end-2), '_{trial\_onset}'];
    end
    ylabel(t);
    if plot_idx == 4
        title('Fisher z-transformed Spearman rank correlation');
    end
end


