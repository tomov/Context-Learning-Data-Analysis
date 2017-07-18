% Compare RDMs from different ROIs with model RDMs
%

%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only

% Get the neural RDMs
%
Neural = rdms_get_neural(data, metadata, which_rows);
showRDMs(Neural, 1);

% Get the model RDMs
%
Model = rdms_get_model(data, metadata, which_rows);
showRDMs(Model, 2);



%% Compute second-order RDM
%
% compare each neural and model RDMs for subject separately using
% Spearman's rank coefficient, then find which ones are significant
%
tic;

% Some knobs / params
%
all_vs_all = false; % #KNOB compare neural vs. neural and model vs. model representations too

control_model_idxs = [8, 12]; % #KNOB control for time and run
assert(isequal(Model(8).name, 'time'));
assert(isequal(Model(12).name, 'run'));

% Compute second-order RDM
%

if all_vs_all
    % all RDMs vs. all RDMs
    %
    % Notice we can't do LME here
    %
    rows = [Neural Model];
    cols = rows;
    [table_Rho, table_H, table_P, all_subject_rhos] = rdms_second_order(metadata, rows, cols, [], false, [], []);
else
    % Neural RDMs vs. Model RDMs
    %
    % It is important that rows = Neural and cols = Model for the visualizations
    % and the linear mixed effects modeling
    %
    %lme_neural_idxs = [1 2 5 11 12 14 15 18 24 25]; % which ROIs to consider for LME
    lme_neural_idxs = 1:numel(Neural);
    lme_model_idxs = [1 3 18 39 41 46]; % which models to consider to LME
    rows = Neural;
    cols = Model;
    [table_Rho, table_H, table_P, all_subject_rhos, lme] = rdms_second_order(metadata, rows, cols, control_model_idxs, true, lme_model_idxs, lme_model_idxs);
end


Rho = array2table(table_Rho, 'RowNames', {rows.name}, 'VariableNames', {cols.name});
H = array2table(table_H, 'RowNames', {rows.name}, 'VariableNames', {cols.name});
P = array2table(table_P, 'RowNames', {rows.name}, 'VariableNames', {cols.name});

toc;





%% Visualize full analysis
% Show full second-order RSA matrix with correlation coefficients + another one with the p-values
%

tabs = {table_Rho, table_P};
titles = {'Representational similarity match (Spearman rank correlation)', 'P-value (one-sample t-test)'};

for i = 1:numel(tabs)
    figure;
    if i == 2
        t = tabs{i};
        if all_vs_all
            imagesc(log10(t), [-16 0.1]);
        else
            imagesc(log10(t));
        end    
        c = colorbar;
        yt = get(c, 'Ytick');
        set(c, 'YTickLabel', arrayfun(@(x) ['10\^', num2str(x)], yt, 'UniformOutput', false));
    else
        if all_vs_all
            imagesc(tabs{i}, [0 0.4]);
        else
            imagesc(tabs{i});
        end
        %cols=colorScale([0 0.5 1; 0.5 0.5 0.5; 1 0 0],256);
        %colormap(cols); colorbar;
        colorbar;
    end

    xticklabels({cols.name});
    set(gca, 'xtick', 1:numel({cols.name}));
    xtickangle(60);
    xlabel('Neural model');

    yticklabels({rows.name});
    set(gca, 'ytick', 1:numel({rows.name}));
    ylabel('ROI_{event}, t = trial onset, f = feedback onset');
    
    title(titles{i});
end













%% Show the significant positive correlations
%
bonferroni = false; % #KNOB -- whether to use Bonferroni correction for the p-values

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
    if linear_mixed_effects
        % LMEs
        %
        col_idxs = [];
        means = [];
        sems = [];
        ps = [];
        for col_idx = 1:numel(cols)
            if isempty(lmes{row_idx, col_idx}), continue; end
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
    if linear_mixed_effects
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
        text(xs(i) - length(stars) * 0.3 / 2, means(i) + sems(i) * 1.2, stars, 'Rotation', 45);
    end

    % Put the ROI / model names and figure title
    %
    set(gca, 'xtick', xs);
    models = {Model(col_idxs).name};
    labels = {};
    for j = 1:numel(col_idxs)
        labels{j} = sprintf('p < 10^{%d}', round(log10(ps(j))));
    end
    xticklabels(models);
    xtickangle(55);

    t = Neural(row_idx).name;
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













%% Interesting visualizations
% show select neural / model pairs that I think are interesting, like bar
% plots
%
assert(~all_vs_all); % only works for neural vs. model comparisons

neural_idxs = [1, 5, 2, 3, 4]; % ROIs to plot at trial onset
neural_idxs = [neural_idxs; neural_idxs + numel(Neural)/2]; % at feedback too
neural_idxs = neural_idxs(:)';

model_idxs = [1:7, 18, 20, 26, 29, 30, 31, 32, 15, 35, 8];

relative_to_time = false; % #KNOB whether to do the analysis/plots relative to time
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
    models = {Model(model_idxs).name};
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
    yticklabels({Neural(neural_idx).name});
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












%% Chan et al.-like bar plots for their models
%
assert(~all_vs_all); % only works for neural vs. model comparisons

neural_idxs = [1:(numel(Neural) / 2)]; % ROIs to plot -- trial_onset
neural_idxs = [1:(numel(Neural) / 2)] + (numel(Neural) / 2); % ROIs to plot -- feedbac_onset
model_idxs = [1:7, 28, 9, 10, 11, 8]; % models to plot -- Chan et al.
relative_to_time = false; % #KNOB whether to do the analysis/plots relative to time
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
    rois = {Neural(neural_idxs).name};
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
    yticklabels({Model(model_idx).name});
    %ytickangle(30);
    if i == 1
        title('Fisher z-transformed Spearman rank correlation');
    end
end



