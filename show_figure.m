function show_figure(figure_name)

sem = @(x) std(x) / sqrt(length(x));

% Generate figures from the paper.
%
% figure_name = which figure to show, e.g. Figure_3A
%

rng default; % be consistent please

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
    
    case 'Figure_curves'
        figure('pos', [100 100 693+20 320] * 3/4);
        
        %
        % Learning curves
        % 
        %figure;
            
        for group = 1:2
            if group == 1
                % pilot subjects
                [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
            else
                % fmri subjects
                [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
            end
            simulated = simulate_subjects(data, metadata, params, which_structures);        
            
            % average learning curve
            %
            human_correct = [];
            model_correct = [];
            
            human_correct_sem = [];
            model_correct_sem = [];

            for n = 1:metadata.trainingTrialsPerRun

                human_corr_n = []; % accuracy on trial n for each subject, averaged across blocks
                model_corr_n = [];
                for who = metadata.subjects
                    which = data.which_rows & data.isTrain & data.trialId == n & strcmp(data.participant, who);
                    assert(sum(which) == 9);
                    
                    subj_corr_n = strcmp(data.response.keys(which), data.corrAns(which)); % accuracy on trial n for subject who, averaged across blocks
                    sim_subj_corr_n = strcmp(simulated.keys(which), data.corrAns(which)); % accuracy on trial n for model for subject who, averaged across blocks
                    
                    human_corr_n = [human_corr_n, mean(subj_corr_n)];
                    model_corr_n = [model_corr_n, mean(sim_subj_corr_n)];
                end                
                
                human_correct = [human_correct mean(human_corr_n)];
                model_correct = [model_correct mean(model_corr_n)];
                
                human_correct_sem = [human_correct_sem sem(human_corr_n)];
                model_correct_sem = [model_correct_sem sem(model_corr_n)];
            end

            subplot(1,2,group);
            plot(model_correct, '.-', 'LineWidth', 2); % == mean(human_correct_all_runs)
            %errorbar(1:metadata.trainingTrialsPerRun, model_correct, model_correct_sem, 'o-', 'LineWidth', 2);
            hold on;
            plot(human_correct, '.-', 'LineWidth', 2); % == mean(model_correct_all_runs)
            %errorbar(1:metadata.trainingTrialsPerRun, human_correct, human_correct_sem, 'o-', 'LineWidth', 2);
            hold off;
            legend({'model', 'subjects'}, 'Position', [0.30 0.05 1 1]);
            %title('Average per-trial accuracy');
            if group == 1
                title('Behavioral pilot');
                text(-4, 1.05, 'A', 'FontSize', 20, 'FontWeight', 'bold')
            else
                title('fMRI');
                text(-4, 1.05, 'B', 'FontSize', 20, 'FontWeight', 'bold')
            end
            xlabel('trial #');
            ylabel('accuracy');
            set(gca,'fontsize',13);
        end
        
    
    % statistics regarding learning and performance on the training trials
    %
    case 'Learning_stats'
        
        for group = 1:2
            if group == 1
                % pilot subjects
                [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
                disp('\n\n\n --------------------------------------- PILOT -----------------\n\n');
            else
                % fmri subjects
                [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
                disp('\n\n\n --------------------------------------- fMRI -----------------\n\n');
            end
            
            second_half_corr = [];
            for who = metadata.subjects
                % accuracy during the second half of training, collapsed across
                % blocks and trials
                %
                which = data.which_rows & data.isTrain & data.trialId > (metadata.trainingTrialsPerRun / 2) & strcmp(data.participant, who);
                assert(sum(which) == 90);

                subj_2nd_half_corr = strcmp(data.response.keys(which), data.corrAns(which)); % accuracy on trial n for subject who, averaged across blocks

                second_half_corr = [second_half_corr, mean(subj_2nd_half_corr)];
            end

            disp(second_half_corr);
            fprintf('accuracy during 2nd half of training: %.6f +- %.6f\n', mean(second_half_corr), sem(second_half_corr));
            disp('t-test: is accuracy = 50%');
            [h, p, ci, stats] = ttest(second_half_corr, 0.5)
        end
        
        
    % statistics showing the generalization pattern on the test trials is real
    %
    case 'generalization_stats'
        human_choices = [];      
        
        % irrelevant condition
        % does old cue cause sickness more than the new cue, in either
        % context?
        %
        condition = 'irrelevant';
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
        
        old_cue = [x1c1; x1c3];
        new_cue = [x3c1; x3c3];
        
        old_cue = mean(old_cue); % collapse across blocks & context (c1 or c3)
        new_cue = mean(new_cue);
        assert(numel(old_cue) == metadata.N);
        assert(numel(new_cue) == metadata.N);

        disp('IRRELEVANT condition t-test: is old cue more predictive of sickness than new cue?');
        [h, p, ci, stats] = ttest2(old_cue, new_cue)
        
        % modulatory condition
        % does old cue-context pair cause sickness more than the other
        % three pairs?
        %
        condition = 'modulatory';
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
        
        old_pair = [x1c1];
        new_pairs = [x3c1; x1c3; x3c3];
        
        old_pair = mean(old_pair); % collapse across blocks & pair (c1 or c3)
        new_pairs = mean(new_pairs);
        assert(numel(old_pair) == metadata.N);
        assert(numel(new_pairs) == metadata.N);

        disp('MODULATORY condition t-test: is old pair more predictive of sickness than new pairs?');
        [h, p, ci, stats] = ttest2(old_pair, new_pairs)        
        

        % additive condition
        % does old context cause sickness more than the new context, for either
        % cues?
        %
        condition = 'additive';
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
        
        old_context = [x1c1; x3c1];
        new_context = [x1c3; x3c3];
        
        old_context = mean(old_context); % collapse across blocks & context (c1 or c3)
        new_context = mean(new_context);
        assert(numel(old_context) == metadata.N);
        assert(numel(new_context) == metadata.N);

        disp('ADDITIVE condition t-test: is old context more predictive of sickness than new context?');
        [h, p, ci, stats] = ttest2(old_context, new_context)
        
        
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
        
        text(0.1, 1.25, 'A', 'FontSize', 20, 'FontWeight', 'bold')

        
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
        
        % Plot model choices probabilities
        %
        h = bar(model_means);

        hold on;
        
        % plot both pilot and fmri subjects' choices
        %
        group_x_offs = [-0.02, 0.02];
        group_color = {[0.5 0.5 0.5], [0 0 0]};
        for group = 1:2
            if group == 1
                % pilot subjects
                [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
            else
                % fmri subjects
                [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
            end
            simulated = simulate_subjects(data, metadata, params, which_structures);        
            
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
            errorbar(xs + group_x_offs(group), human_means, human_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', group_color{group}, 'LineWidth', 1, 'Color', group_color{group}, 'AlignVertexCenters', 'off');
            % the markers in the errorbar plot are misaligned from the error
            % bars -- this hack adjusts them
            xs_adjust = xs + 0.003;
            %{
            xs_adjust(1) = xs_adjust(1) + 0.0035;
            xs_adjust(2) = xs_adjust(2) + 0.003;
            xs_adjust(3) = xs_adjust(3) + 0.003;
            xs_adjust(4) = xs_adjust(4) + 0.0025;
            xs_adjust(5) = xs_adjust(5) + 0.004;
            xs_adjust(6) = xs_adjust(6) + 0.004;
            xs_adjust(7) = xs_adjust(7) + 0.0025;
            xs_adjust(8) = xs_adjust(8) + 0.0025;
            xs_adjust(9) = xs_adjust(9) + 0.003;
            xs_adjust(10) = xs_adjust(10) + 0.004;
            xs_adjust(11) = xs_adjust(11) + 0.003;
            xs_adjust(12) = xs_adjust(12) + 0.0025;
            %}
            plot(xs_adjust + group_x_offs(group), human_means, 'o', 'MarkerSize', 5, 'MarkerFaceColor', group_color{group}, 'Color', group_color{group});
        end
        
        
        
        hold off;
        
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Choice probability');
        legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'}, 'Position', [0.07 -0.095 1 1]);
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        
        %print(gcf, 'untitled.pdf', '-dpdf', '-bestfit');
        
        text(0.1, 1.25, 'B', 'FontSize', 20, 'FontWeight', 'bold')

        
        
    case 'Figure_3_stats'        
        
        %
        % Correlate model with subject choices
        %
        load(fullfile('results', 'fit_params_results.mat'), 'results', 'results_options');        
        
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
        bic = results(1).bic; % TODO FIXME the params are close but not quite there
        fprintf('Total log likelihood = %f; BIC = %f (on whatever it was fit)\n', total_loglik, bic);

        %
        % Correlate single-causal-structure models with human data
        %
        
        results_idxs = [2 3 4];
        which_structuress = {[1 0 0 0], [0 1 0 0], [0 0 1 0]};
        structure_names = {'M_1', 'M_2', 'M_3'};
        
        for i = 1:3
            params = results(results_idxs(i)).x;
            options = results_options(results_idxs(i));
            bic = results(results_idxs(i)).bic;
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
            fprintf('Total log likelihood for %s = %f; BIC = %f (on whatever it was fit)\n', structure_names{i}, total_loglik, bic);
        end
        
        assert(false, 'THESE ARE DEPRECATED -- use the stats from plot_behavior.m; you have to change which_structures e.g. to [1 0 0 0] for M1 only');
        
        
        
    % same as Figure 3 except for pilot data
    %
    case 'Figure_S1'
        
        [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
        simulated = simulate_subjects(data, metadata, params, which_structures);
        
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
        
        text(0.1, 1.25, 'A', 'FontSize', 20, 'FontWeight', 'bold')
        
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
        errorbar(xs, human_means, human_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        % the markers in the errorbar plot are misaligned from the error
        % bars -- this hack adjusts them
        xs_adjust = xs;
        xs_adjust(1) = xs_adjust(1) + 0.0035;
        xs_adjust(2) = xs_adjust(2) + 0.003;
        xs_adjust(3) = xs_adjust(3) + 0.003;
        xs_adjust(4) = xs_adjust(4) + 0.0025;
        xs_adjust(5) = xs_adjust(5) + 0.004;
        xs_adjust(6) = xs_adjust(6) + 0.004;
        xs_adjust(7) = xs_adjust(7) + 0.0025;
        xs_adjust(8) = xs_adjust(8) + 0.0025;
        xs_adjust(9) = xs_adjust(9) + 0.003;
        xs_adjust(10) = xs_adjust(10) + 0.004;
        xs_adjust(11) = xs_adjust(11) + 0.003;
        xs_adjust(12) = xs_adjust(12) + 0.0025;
        plot(xs_adjust, human_means, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0 0 0], 'Color', [0 0 0]);
        hold off;
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('Choice probability');
        legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'}, 'Position', [0.07 -0.095 1 1]);
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        
        
        text(0.1, 1.25, 'B', 'FontSize', 20, 'FontWeight', 'bold')

        
        %print(gcf, 'untitled.pdf', '-dpdf', '-bestfit');
        
        
        
    case 'searchlight_posterior'
        bspmview('rdms/betas_smooth/searchlight_tmap_posterior_feedback_onset.nii', '../neural/mean.nii');
        
    case 'searchlight_prior'
        bspmview('rdms/betas_smooth/searchlight_tmap_prior_trial_onset.nii', '../neural/mean.nii');
        
    case 'KL_structures'
        ccnl_view(context_expt(), 154, 'KL_structures');
        
    case 'KL_weights'
        ccnl_view(context_expt(), 154, 'KL_weights');
        
    case 'KL_structures - KL_weights'
        ccnl_view(context_expt(), 154, 'KL_structures - KL_weights');

    case 'KL_weights - KL_structures'
        ccnl_view(context_expt(), 154, 'KL_weights - KL_structures');
        

    case 'glm154'
        figure('pos', [100 100 693+20 492]);
        
        fontsize = 12;
        
        %
        % Top left: structures 
        %

        subplot(2, 2, 1);
        
        imshow('images/KL_structures_pos.png'); % from GLM 154
        title('Causal structure update', 'FontSize', fontsize);
       
        %
        % Top right: weights 
        %
        
        subplot(2, 2, 2);
        
        imshow('images/KL_weights_pos.png'); % from GLM 154
        title('Associative weights update', 'FontSize', fontsize);
        
        %
        % Bottom left: contrast
        %
        
        subplot(2, 2, 3);
        
        imshow('images/KL_structures-KL_weights.png'); % from GLM 154
        title({'Causal structure update >'; 'associative weights update'}, 'FontSize', fontsize);
        
        %
        % Bottom right: KL structures ~ test choice log likelihood
        %

        subplot(2, 2, 4);
        
        load results/betas_to_behavior_glm154_KL_structures_sphere_KL_structures.mat
        %which = 1:numel(region);
        assert(numel(region) == 14);
        which = [1:7 12 14]; % exclude motor and visual areas
        % exclude motor & cerebellum, also reorder them
       % which = [1 4 2 6 5 3];
        
        region = region(which);
        means = means(which);
        sems = sems(which);
        ps = ps(which);
        mni = mni(which, :);

        h = bar(means, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
        xs = h(1).XData;
        hold on;
        errorbar(xs, means, sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        for i = 1:numel(xs)
            if ps(i) < 0.05
                assert(ps(i) >= 0.01);
                text(xs(i) - 0.08, means(i) + sems(i) * 1.1, '*', 'FontSize', fontsize);
            end
        end
        hold off;

        set(gca, 'xtick', xs);
        ylabel('Fisher z-transformed r');
        
        neural_names = {};
        for i = 1:numel(region)
            hemisphere = [];
            if mni(i, 1) < 0
                hemisphere = 'L';
            else
                hemisphere = 'R';
            end
            roi = aal2_label_to_roi_name(region{i});
            comma = strfind(roi, ','); % trim the super long names
            if ~isempty(comma)
                roi = roi(1:comma(1)-1);
            end
            neural_names{i} = [roi, ' (', hemisphere, ')'];
        end
        % screw that, just hardcode the names
        neural_names = { ...
            'IFG pars orbitalis (R)', ...
            'Angular gyrus (R)', ...
            'IFG pars opercularis (R)', ...
            'Frontal Pole (L)', ...
            'Angular gyrus (L)', ...
            'Middle orbital gyrus (R)', ...
            'IFG pars opercularis (L)', ...
            'Inferior temporal gyrus (R)', ...
            'Lingual gyrus (R)', ...
            'Middle temporal gyrus (L)'};
        
        xticklabels(neural_names);
        xtickangle(30);
        title('Structure updating predicts test choices', 'FontSize', fontsize);

        % label subplots A,B,C,D
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        axes(ax1);
        text(0.1, 0.95, 'A', 'FontSize', 20, 'FontWeight', 'bold');
        text(0.51, 0.95, 'B', 'FontSize', 20, 'FontWeight', 'bold');
        text(0.1, 0.47, 'C', 'FontSize', 20, 'FontWeight', 'bold');
        text(0.51, 0.47, 'D', 'FontSize', 20, 'FontWeight', 'bold');
        

    case 'searchlight'
        
        figure('pos', [100 100 693+20 492]);
        
        fontsize = 12;
        
        %
        % Top left: structures 
        %

        subplot(1, 2, 1);
        
        imshow('images/searchlight_prior.png');
        title('Causal structure prior', 'FontSize', fontsize);
        
        g = subplot(1, 2, 2);
        
        imshow('images/searchlight_posterior.png');
        title('Causal structure posterior', 'FontSize', fontsize);
        p = get(g, 'position');
        p([3 4]) = p([3 4]) * 1.025;
        p(2) = p(2) - 0.01;
        set(g, 'position', p);
        
        % label subplots A,B,C,D
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        axes(ax1);
        text(0.1, 0.7, 'A', 'FontSize', 20, 'FontWeight', 'bold');
        text(0.54, 0.7, 'B', 'FontSize', 20, 'FontWeight', 'bold');

        
        
        







        
    case 'Figure_4A_Neuron'
        ccnl_view(context_expt(), 123, 'surprise');
        
    case 'Figure_4C_Neuron'
        ccnl_view(context_expt(), 148, 'KL_weights');
        

    case 'Figure_4A_stats_Neuron'        
        
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
                
        
        
        
        
    case 'Figure_4_Neuron'
        figure;
        
        %
        % Top row: weights
        %

        % Plot saved brain activation map for weights KL
        %
        subplot(2, 2, 3);
        
        %imshow('images/kl-weights.png'); % from glm 148, 'KL_weights'
        imshow('images/kl-weights-surface.png'); % from glm 148, 'KL_weights'
        title('Associative weights', 'FontSize', 10);
        
        % Plot fisher Z-transformed correlation coefficients between
        % per-run log likelihood of subject test choices and temporal beta for
        % weights KL
        %
        load(fullfile('results', 'kl_weights.mat')); % computed by kl_weights.m

        subplot(3, 4, 11);
        
        assert(isequal(rois{1}, 'Temporal_Inf_L'));
        rs1 = all_fisher_rs{1};
        [rs1, subj_ids] = sort(rs1);
        
        bar(rs1, 'EdgeColor', 'none');
        set(gca, 'XTick', 1:1:20);
        %xticklabels(subj_ids);
        xticklabels({[]});
        
        xlabel('Subject', 'FontSize', 10);
        ylabel('Correlation with behavior', 'FontSize', 10);
        xlim([0 21]);
        ylim([-0.75 1.1]);
        title('L Inf Temporal Gyrus', 'FontSize', 10);
        set(gca, 'xtick', []);

        % same for parietal betas
        %
        subplot(3, 4, 12);
        
        assert(isequal(rois{2}, 'Parietal_Inf_L'));
        rs2 = all_fisher_rs{2};
        [rs2, subj_ids] = sort(rs2);
        
        bar(rs2, 'EdgeColor', 'none');
        set(gca, 'XTick', 1:1:20);
        %xticklabels(subj_ids);
        xticklabels({[]});
        
        xlabel('Subject', 'FontSize', 10);
        %ylabel('Correlation with behavior');
        xlim([0 21]);
        ylim([-0.75 1.1]);
        title('L Angular Gyrus', 'FontSize', 10);
        set(gca, 'xtick', []);
        
        %
        % Bottom row: posterior
        %
        
        % Plot saved brain activation map for posterior KL
        %
        subplot(2, 2, 1);
        
        %imshow('images/kl-posterior.png'); % from glm 123, 'surprise'
        imshow('images/kl-structures-surface.png'); % from glm 123, 'surprise'
        title('Causal structures', 'FontSize', 10);
        
        % Plot fisher Z-transformed correlation coefficients between
        % per-run log likelihood of subject test choices and AG beta for
        % posterior KL
        %
        load(fullfile('results', 'kl_posterior.mat')); % computed by kl_weights.m

        subplot(2, 2, 2);
        
        assert(isequal(rois{1}, 'Parietal_Sup_R'));
        rs = all_fisher_rs{1};
        [rs, subj_ids] = sort(rs);
        
        bar(rs, 'EdgeColor', 'none');
        set(gca, 'XTick', 1:1:20);
        %xticklabels(subj_ids);
        xticklabels({[]});
        
        xlabel('Subject', 'FontSize', 10);
        ylabel('Correlation with behavior', 'FontSize', 10);
        xlim([0 21]);
        ylim([-0.75 1.1]);
        title('R Angular Gyrus', 'FontSize', 10);
        set(gca, 'xtick', []);
        %set(gca,'fontsize',13);
        
        %print(gcf, 'Figure_4B.png', '-dpng', '-r300');
        %print(gcf, 'images/fmri-results.pdf', '-dpdf', '-bestfit');
        
        
        
        
        
        
        
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
