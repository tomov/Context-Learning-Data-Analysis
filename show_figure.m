function show_figure(figure_name)

% Generate figures from the paper.
%
% figure_name = which figure to show, e.g. Figure_3A
%

% set nicer colors
C = linspecer(3);
set(0,'DefaultAxesColorOrder',C)

filename = 'show_figure.mat';

% Try to load cached figure data
%
if exist(filename, 'file') ~= 2
    % Load data
    %
    [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

    % Load parameters
    %
    load(fullfile('results', 'fit_params_results_fmri_random_effects_20_nstarts_5_prior.mat'), 'results', 'results_options');
    params = results(1).x;
    options = results_options(1);
    
    % OVERRIDE -- use pilot params as before
    %
    params = [0.1249 2.0064];
    options.isFmriData = false;
    options.fixedEffects = true;
    
    disp('Using parameters:');
    disp(params);
    disp('generated with options:');
    disp(options);
    % safeguards
    %assert(options.isFmriData == true);
    %assert(~options.fixedEffects);
    assert(isequal(options.which_structures, [1 1 1 0]));
    which_structures = logical([1 1 1 0]);

    % Run the model with the parameters
    %
    simulated = simulate_subjects(data, metadata, params, which_structures);
    
    fprintf('Saving data to %s\n', filename);
    save(filename, 'data', 'metadata', 'simulated', 'results', 'params', 'options', 'which_structures'); % cache to file
else
    fprintf('Loading data from %s\n', filename);
    load(filename, 'data', 'metadata', 'simulated', 'results', 'params', 'options', 'which_structures');
end

% Plot figure(s)
%
switch figure_name
    
    case 'Figure_3'
        
        %
        % Figure 3A: Posterior probabilities over structures in each condition
        %
        
        figure;
        %set(handle, 'Position', [500, 500, 450, 200])
        
        subplot(2, 1, 1);

        P_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.isTrain & data.trialId == 20 & strcmp(data.contextRole, condition);
            
            P = simulated.P(which_rows, which_structures);             
            P_means = [P_means; mean(P, 1)];
        end

        bar(P_means);
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Posterior probability');
        legend({'M1', 'M2', 'M3'}, 'Position', [0.15 0.3 1 1]);
        ylim([0 1.1]);
        set(gca,'fontsize',13);    
        
        %
        % Figure 3B: Choice probabilities on test trials for model vs. humans
        %
        
        subplot(2, 1, 2);
        
        % Choice probabilities for model
        %
        model_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & ~data.isTrain & strcmp(data.contextRole, condition);
            
            x1c1 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 0);
            x1c3 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 2);
            x3c1 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 0);
            x3c3 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 2);

            model_means = [model_means; mean(x1c1) mean(x1c3) mean(x3c1) mean(x3c3)];
        end

        % Choice probabilities for human subjects
        %
        human_choices = [];        
        for condition = metadata.contextRoles
            which_rows = data.which_rows & ~data.isTrain & strcmp(data.contextRole, condition);
            
            x1c1 = []; x1c3 = []; x3c1 = []; x3c3 = [];
            for who = metadata.subjects
                x1c1 = [x1c1, data.chose_sick(which_rows & data.cueId == 0 & data.contextId == 0 & strcmp(data.participant, who))];
                x1c3 = [x1c3, data.chose_sick(which_rows & data.cueId == 0 & data.contextId == 2 & strcmp(data.participant, who))];
                x3c1 = [x3c1, data.chose_sick(which_rows & data.cueId == 2 & data.contextId == 0 & strcmp(data.participant, who))];
                x3c3 = [x3c3, data.chose_sick(which_rows & data.cueId == 2 & data.contextId == 2 & strcmp(data.participant, who))];
            end
            assert(isequal(size(x1c1), [metadata.runsPerContext metadata.N]));
            assert(isequal(size(x1c3), [metadata.runsPerContext metadata.N]));
            assert(isequal(size(x3c1), [metadata.runsPerContext metadata.N]));
            assert(isequal(size(x3c3), [metadata.runsPerContext metadata.N]));
            
            human_choices = [human_choices; mean(x1c1); mean(x1c3); mean(x3c1); mean(x3c3)];
        end
        
        [human_sems, human_means] = wse(human_choices');
        
        % Plot model choices probabilities
        %
        h = bar(model_means);
        
        % Get x-coordinates of bar centers and "flatten" the human choice
        % data so we can overlay it nicely
        %
        xs = sort([h(1).XData + h(1).XOffset, ...
                   h(2).XData + h(2).XOffset, ...
                   h(3).XData + h(3).XOffset, ...
                   h(4).XData + h(4).XOffset]);
        human_means = human_means'; human_means = human_means(:);
        human_sems = human_sems'; human_sems = human_sems(:);        
        
        % Plot human choices
        %
        hold on;
        errorbar(xs, human_means, human_sems, 'o', 'LineWidth', 2, 'Color', [0 0 0]);
        hold off;
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Choice probability');
        legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'}, 'Position', [0.07 -0.095 1 1]);
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        
        
    case 'Figure_3_stats'        
        
        %
        % Correlate model with subject choices
        %
        
        s = simulated.pred(data.which_rows);
        d = data.chose_sick(data.which_rows);
        [r, p] = corrcoef([s, d]);
        r = r(1,2);
        p = p(1,2);
        fprintf('Correlation between model and humans (ALL TRIALS): r = %f, p = %f\n', r, p);
        
        s = simulated.pred(data.which_rows & ~data.isTrain);
        d = data.chose_sick(data.which_rows & ~data.isTrain);
        [r, p] = corrcoef([s, d]);
        r = r(1,2);
        p = p(1,2);
        fprintf('Correlation between model and humans (TEST TRIALS): r = %f, p = %f\n', r, p);
        
        total_loglik = model_likfun(data, metadata, params, which_structures, data.which_rows, false);
        fprintf('Total log likelihood = %f\n', total_loglik);

        %
        % Correlate single-causal-structure models with human data
        %
        
        results_idxs = [2 3 4];
        which_structuress = {[1 0 0 0], [0 1 0 0], [0 0 1 0]};
        structure_names = {'M_1', 'M_2', 'M_3'};
        
        for i = 1:3
            load(fullfile('results', 'fit_params_results.mat'), 'results', 'results_options');
            params = results(results_idxs(i)).x;
            options = results_options(results_idxs(i));
            %assert(options.isFmriData == true);
            %assert(~options.fixedEffects);
            assert(isequal(options.which_structures, which_structuress{i}));

            simulated = simulate_subjects(data, metadata, params, which_structuress{i});

            s = simulated.pred(data.which_rows);
            d = data.chose_sick(data.which_rows);
            [r, p] = corrcoef([s, d]);
            r = r(1,2);
            p = p(1,2);

            s_test = simulated.pred(data.which_rows & ~data.isTrain);
            d_test = data.chose_sick(data.which_rows & ~data.isTrain);
            [r_test, p_test] = corrcoef([s_test, d_test]);
            r_test = r_test(1,2);
            p_test = p_test(1,2);
            
            fprintf('$r = %.4f$ for all trials ($r = %.4f$ for the test trials) for $%s$\n', r, r_test, structure_names{i});            

            total_loglik = model_likfun(data, metadata, params, which_structuress{i}, data.which_rows, false);
            fprintf('Total log likelihood for %s = %f\n', structure_names{i}, total_loglik);
        end
        
        
        
    case 'Figure_4A'
        % Slices showing KL divergence activations [use the model with the error regressor]
        %
        ccnl_view(context_expt(), 123, 'surprise - wrong');
        

    case 'Figure_4A_stats'        
        
        %
        % Correlate D_KL and wrong answers
        %
        s = simulated.surprise(data.which_rows & data.isTrain);
        model_keys = simulated.keys(data.which_rows & data.isTrain);
        corr_ans = data.corrAns(data.which_rows & data.isTrain);
        model_corr = strcmp(model_keys, corr_ans);
        d = ~model_corr;
        [r, p] = corrcoef([s, d]);
        r = r(1,2);
        p = p(1,2);
        fprintf('Correlation D_KL and model errors (TRAINING TRIALS): r = %f, p = %f\n', r, p);
        
        s = simulated.surprise(data.which_rows & data.isTrain);
        d = ~data.response.corr(data.which_rows & data.isTrain);
        [r, p] = corrcoef([s, d]);
        r = r(1,2);
        p = p(1,2);
        fprintf('Correlation D_KL and human errors (TRAINING TRIALS): r = %f, p = %f\n', r, p);
                
        
    case 'Figure_4B'
        % Plot fisher Z-transformed correlation coefficients between
        % per-run log likelihood of subject test choices and AG beta
        %
        
        load('kl_analysis.mat'); % as output by kl_structure_learning.m

        figure;
        
        assert(isequal(rois{1}, 'Angular_R'));
        rs = fisher_all_rs(1,:);
        [rs, subj_ids] = sort(rs);
        
        bar(rs);
        set(gca, 'XTick', 1:1:20);
        xticklabels(subj_ids);
        xlabel('Subject #');
        ylabel('Fisher z-transformed correlation coefficient');
        xlim([0 21]);
        set(gca,'fontsize',13);
        
        
        
    case 'Figure_4B_OLD'     
        
        % Plot showing the least-squares lines relating AG KL beta to the structure learning effect, one line per subject, plus a thicker line showing the average linear relationship.    
        %
        z_scored = false;
        
        load('kl_analysis.mat'); % as output by kl_structure_learning.m

        handle = figure;
        set(handle, 'Position', [500, 500, 200, 200])
        
        assert(strcmp(rois{1}, 'Angular_R'));
        KL_betas_AG = kl_betas(:, :, 1);
        slopes = nan(1, n_subjects);
        intercepts = intercepts(1, n_subjects);
        xs = nan(n_runs, n_subjects);
        
        for subj_idx = 1:n_subjects
            %x = structure_learnings(:, subj_idx);   <-- not good with timeouts
            x = test_liks(:, subj_idx);
            y = KL_betas_AG(:, subj_idx);

            if z_scored
                x = zscore(x);
                y = zscore(y);
            end
            xs(:, subj_idx) = x;
            
            fit = polyfit(x, y, 1);
            slopes(subj_idx) = fit(1);
            intercepts(subj_idx) = fit(2);
            
            hold on;
            yfit = polyval(fit, x);
            plot(x, yfit, 'Color', [0.5 0.5 0.5]);
            hold off;
        end
        
        %x = [-2 0];
        %hold on;
        %for subj_idx = 1:n_subjects
        %    plot(x, x * slopes(subj_idx) + intercepts(subj_idx), 'Color', 'blue');
        %end
        
        max_x = max(xs(:));
        
        hold on;
        x_limits = [-1.6 max_x];
        plot(x_limits, x_limits * mean(slopes) + mean(intercepts), 'LineWidth', 2, 'Color', 'blue');
        hold off;
        
        xlim(x_limits);
        title('Right Angular Gyrus');
        if z_scored
            xlabel('Test choice likelihood (Z-scored)');
            ylabel('Beta for KL divergence (Z-scored)');
            ylim([-1 1]);
        else
            xlabel('Test choice likelihood');
            ylabel('Beta for KL divergence');
            ylim([-30 35]);
        end

end


%{
Nice, I think we can build a story around this:
(1) People do structure learning, as demonstrated by our behavior and modeling.
(2) Neural correlates of posterior updating in several regions, including angular gyrus and lateral OFC.
(3) Variability in Bayesian updating within the AG predicts variability in performance at test.

I think we should write this up for a short format journal, starting with Current Biology. I will write a presubmission inquiry to the editor. The next thing we should be thinking about is arranging the figures and then writing the paper around them. Can you organize all the relevant data on dropbox?

In terms of figures, here's what I'm thinking:

Figure 1
Hypothesis space of causal structures (Figure 1 from my PBR paper).

Figure 2
Experimental design: sequence of events in a trial

Figure 3
Model predictions and behavioral results
(A) Posterior probabilities over structures in each condition
(B) Test performance: model predictions as bars, human data as dots + standard error bars overlaid on the bars

Figure 4
Neural results
(A) Slices showing KL divergence activations [use the model with the error regressor]
(B) Plot showing the least-squares lines relating AG KL beta to the structure learning effect, one line per subject, plus a thicker line showing the average linear relationship.
%}
