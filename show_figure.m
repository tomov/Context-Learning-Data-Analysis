function show_figure(figure_name)

utils; % include some nifty lambdas
sem = @(x) std(x) / sqrt(length(x));

% Generate figures from the paper.
%
% figure_name = which figure to show, e.g. Figure_3A
%

rng default; % reproducibility

% set nicer colors
C = linspecer(3);
set(0,'DefaultAxesColorOrder',C);
   


% Plot figure(s)
%
switch figure_name
    
    case 'fig:curves'
        %
        % Learning curves, model vs. subjects
        % 
        figure('pos', [100 100 693+20 320] * 3/4);
     
        % pilot 
        %[data, metadata, simulated] = simulate_subjects_helper(false);
        [data, metadata, simulated] = simulate_subjects_helper(false, fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat'), 1, [1 1 0 1 0]);
        subplot(1,2,1);
        plot_curves_helper(data, metadata, simulated);
        title('Behavioral pilot');
        text(-4, 1.05, 'A', 'FontSize', 20, 'FontWeight', 'bold')

        % fmri
        %[data, metadata, simulated] = simulate_subjects_helper(true);
        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat'), 1, [1 1 0 1 0]);
        subplot(1,2,2);
        plot_curves_helper(data, metadata, simulated);
        title('fMRI');
        text(-4, 1.05, 'B', 'FontSize', 20, 'FontWeight', 'bold')

    case 'fig:curves_by_condition'
        %
        % Learning curves, model vs. subjects
        % 
        figure('pos', [100 100 693+20 320] * 3/4);
     
        % pilot 
        [data, metadata, simulated] = simulate_subjects_helper(false, fullfile('results', 'fit_params_results_M1M2M3_25nstarts.mat'), 1, [1 1 1 0 0]);
        subplot(2,3,1);
        plot_curves_helper(data, metadata, simulated, data.which_rows & strcmp(data.contextRole, 'irrelevant'));
        title('Behavioral pilot: irr');
        subplot(2,3,2);
        plot_curves_helper(data, metadata, simulated, data.which_rows & strcmp(data.contextRole, 'modulatory'));
        title('mod');
        subplot(2,3,3);
        plot_curves_helper(data, metadata, simulated, data.which_rows & strcmp(data.contextRole, 'additive'));
        title('add');
        text(-4, 1.05, 'A', 'FontSize', 20, 'FontWeight', 'bold')

        % fmri
        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_M1M2M3_25nstarts.mat'), 1, [1 1 1 0 0]);
        subplot(2,3,4);
        plot_curves_helper(data, metadata, simulated, data.which_rows & strcmp(data.contextRole, 'irrelevant'));
        title('fMRI: irr');
        subplot(2,3,5);
        plot_curves_helper(data, metadata, simulated, data.which_rows & strcmp(data.contextRole, 'modulatory'));
        title('mod');
        subplot(2,3,6);
        plot_curves_helper(data, metadata, simulated, data.which_rows & strcmp(data.contextRole, 'additive'));
        title('add');
        text(-4, 1.05, 'B', 'FontSize', 20, 'FontWeight', 'bold')
           
    % statistics regarding learning and performance on the training trials
    %
    case 'stats:learning'

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

            if group == 1 
                % pilot subjects
                report.pilot_mean_second_half_corr = mean(second_half_corr) * 100;
                report.pilot_sem_second_half_corr = sem(second_half_corr) * 100;
                report.pilot_t_df = stats.df;
                report.pilot_t_tstat = stats.tstat;
                report.pilot_t_p = p;
                report.pilot_t_p_less_than_pow10 = ceil(log10(p));
            else
                % fmri subjects
                report.fmri_mean_second_half_corr = mean(second_half_corr) * 100;
                report.fmri_sem_second_half_corr = sem(second_half_corr) * 100;
                report.fmri_t_df = stats.df;
                report.fmri_t_tstat = stats.tstat;
                report.fmri_t_p = p;
                report.fmri_t_p_less_than_pow10 = ceil(log10(p));
            end
        end

        paragraph = 'Average accuracy during the second half of training was $%.1f \\pm %.1f\\%%$ ($t_{%d} = %.1f, p < 10^{%d}$, one-sample \\textit{t}-test against 50\\%%) for the pilot subjects, and $%.1f \\pm %.1f\\%%$ ($t_{%d} = %.1f, p < 10^{%d}$, one-sample \\textit{t}-test against 50\\%%) for the scanned subjects, well above chance.\n\n';
        fprintf(paragraph, ...
                report.pilot_mean_second_half_corr, ...
                report.pilot_sem_second_half_corr, ...
                report.pilot_t_df, ...
                report.pilot_t_tstat, ...
                report.pilot_t_p_less_than_pow10, ...
                report.fmri_mean_second_half_corr, ...
                report.fmri_sem_second_half_corr, ...
                report.fmri_t_df, ...
                report.fmri_t_tstat, ...
                report.fmri_t_p_less_than_pow10);
        
       

    % statistics showing the generalization pattern on the test trials is real
    %
    case 'stats:generalization'
        human_choices = [];      

        [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
        
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

        report.irr_t_df = stats.df;
        report.irr_t_tstat = stats.tstat;
        report.irr_t_p = p;
        report.irr_t_p_pow10 = ceil(log10(p));
        
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
        
        report.mod_t_df = stats.df;
        report.mod_t_tstat = stats.tstat;
        report.mod_t_p = p;
        report.mod_t_p_pow10 = ceil(log10(p));
        

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

        report.add_t_df = stats.df;
        report.add_t_tstat = stats.tstat;
        report.add_t_p = p;
        report.add_t_p_pow10 = ceil(log10(p));


        paragraph = 'The new cue $x_3$, on the other hand, was judged to be much less predictive of sickness in either context ($t_{38} = 9.51, p < 10^{-10}$, paired \\textit{t}-test). Conversely, on blocks during which context acted like another cue (Figure~\\ref{fig:behavior}B, additive training), subjects guessed that both cues would cause sickness in the old context $c_1$ (circle for $x_3 c_1$), but not in the new context $c_3$ ($t_{38} = 11.1, p < 10^{-12}$, paired \\textit{t}-test). Generalization in both of these conditions was different from what one would expect if subjects treated each cue-context pair as a unique stimulus independent from the other pairs, which is similar to the generalization pattern on modulatory blocks (Figure~\\ref{fig:behavior}B, modulatory training). On these blocks, subjects judged that the old cue is predictive of sickness in the old context significantly more compared to the remaining cue-context pairs ($t_{38} = 9.01, p < 10^{-10}$, paired \\textit{t}-test)\n';
        fprintf(paragraph, ...
                report.irr_t_df, ...
                report.irr_t_tstat, ...
                report.irr_t_p_pow10, ...
                report.add_t_df, ...
                report.add_t_tstat, ...
                report.add_t_p_pow10, ...
                report.mod_t_df, ...
                report.mod_t_tstat, ...
                report.mod_t_p_pow10);
       


    % how well does the model correlate with subject test phase choices
    %
    case 'stats:model'
        % order is important
        which_structuress = {[1 1 1 0], [1 0 0 0], [0 1 0 0], [0 0 1 0]};

        for i = 1:numel(which_structuress) 
            which_structures = which_structuress{i};

            [r, p] = get_test_choice_correlations([], i, which_structures);

            switch i
                case 1
                    report.full_r = r;
                    report.full_p = p;
                    report.full_p_pow10 = ceil(log10(p));
                case 2
                    report.M1_r = r;
                    report.M1_p = p;
                    report.M1_p_pow10 = ceil(log10(p));
                case 3
                    report.M2_r = r;
                    report.M2_p = p;
                    report.M2_p_pow10 = ceil(log10(p));
                case 4
                    report.M3_r = r;
                    report.M3_p = p;
                    report.M3_p_pow10 = ceil(log10(p));
            end
        end

        fprintf('Correlation between model and subject test behavior: r = %f, p = %.10f (p = 10^%f)\n', r, p, log10(p));

        paragraph = 'Using parameters fit with data from the behavioral pilot version of the study, the model quantitatively accounted for the generalization pattern on the test trials choices of subjects in the fMRI portion of the study (Figure~\\ref{fig:behavior}B; $r = %.2f, p < 10^{%d}$). As expected, the stimulus-outcome contingencies induced the model to infer a different causal structure in each of the three conditions (Figure~\\ref{fig:behavior}A), leading to the distinct response patterns on the simulated test trials. For comparison, we also ran versions of the model using a single causal structure. Theories corresponding to each of these sub-models have been put forward in the literature as explanations of the role of context during learning, however neither of them has been able to provide a comprehensive account of the behavioral findings on its own. Accordingly, model performance was markedly worse when the hypothesis space was restricted to a single causal structure: the correlation coefficients were $r = %.2f$ for the irrelevant context structure ($M_1$; $p = %.2f$), $r = %.2f$ for the modulatory context structure ($M_2$; $p < %.2f$), and $r = %.2f$ for  the additive context structure ($M_3$; $p < %.4f$).\n';
        fprintf(paragraph, ...
                report.full_r, ...
                report.full_p_pow10, ...
                report.M1_r, ...
                report.M1_p, ...
                report.M2_r, ...
                10^report.M2_p_pow10, ...
                report.M3_r, ...
                10^report.M3_p_pow10);


    case 'tab:models'

        load_cached_values = true;
        cached_file = fullfile('results', 'show_figure_tab_models.mat');

        if load_cached_values
            load(cached_file);
        else

            idx = 0;

            idx = idx + 1;
            models(idx).which_structures = [1 1 1 0]; 
            models(idx).name = 'M1, M2, M3';
            models(idx).params_file = fullfile('results', 'fit_params_results.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 0 0 0]; 
            models(idx).name = 'M1';
            models(idx).params_file = fullfile('results', 'fit_params_results.mat');
            models(idx).params_idx = 2;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [0 1 0 0]; 
            models(idx).name = 'M2';
            models(idx).params_file = fullfile('results', 'fit_params_results.mat');
            models(idx).params_idx = 3;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [0 0 1 0]; 
            models(idx).name = 'M3';
            models(idx).params_file = fullfile('results', 'fit_params_results.mat');
            models(idx).params_idx = 4;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;


            idx = idx + 1;
            models(idx).which_structures = 'simple_Q'; 
            models(idx).name = 'Q learning';
            models(idx).params_file = fullfile('results', 'fit_params_results_simple_q.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\beta = %.2f';
            models(idx).do_include = false;


            idx = idx + 1;
            models(idx).which_structures = 'Q_learning'; 
            models(idx).name = 'Q learning 2';
            models(idx).params_file = fullfile('results', 'fit_params_results_q_learning.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\alpha = %.2f, \\beta = %.2f';
            models(idx).do_include = false;


            idx = idx + 1;
            models(idx).which_structures = [1 1 0 1 0]; 
            models(idx).name = 'M1, M2, M1''';
            models(idx).params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 0 1 0 1]; 
            models(idx).name = 'M1, M2'', M3';
            models(idx).params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
            models(idx).params_idx = 2;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 0 0 1 1]; 
            models(idx).name = 'M1, M2'', M1''';
            models(idx).params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
            models(idx).params_idx = 3;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [0 0 0 1 0]; 
            models(idx).name = 'M1''';
            models(idx).params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
            models(idx).params_idx = 4;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [0 0 0 0 1]; 
            models(idx).name = 'M2''';
            models(idx).params_file = fullfile('results', 'fit_params_results_reviewer2.mat');
            models(idx).params_idx = 5;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 1 1 1 0]; 
            models(idx).name = 'M1, M2, M3, M1''';
            models(idx).params_file = fullfile('results', 'fit_params_results_reviewer2_four.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = 'simple_collins'; 
            models(idx).name = 'Collins 2016, $1<\alpha<2$';
            models(idx).params_file = fullfile('results', 'fit_params_results_simple_collins_5nstarts.mat'); % concentration parameter between 1 and 2
            models(idx).params_idx = 1;
            models(idx).params_format = '\\eta = %.2f, \\beta = %.2f, \\alpha = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = 'flat_collins'; 
            models(idx).name = 'Flat RL';
            models(idx).params_file = fullfile('results', 'fit_params_results_flat_collins_5nstarts.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\eta = %.2f, \\beta = %.2f';
            models(idx).do_include = false;


            % these all have 25 nstarts


            idx = idx + 1;
            models(idx).which_structures = 'simple_collins'; 
            models(idx).name = 'Collins 2016';
            models(idx).params_file = fullfile('results', 'fit_params_results_simple_collins_25nstarts_0-10alpha.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\eta = %.2f, \\beta = %.2f, \\alpha = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = 'simple_collins'; 
            models(idx).name = 'RL + clustering'; % Collins et al 2013, 2016
            models(idx).params_file = fullfile('results', 'fit_params_results_simple_collins_25nstarts_0-10alpha_Q0.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\eta = %.2f, \\beta = %.2f, \\alpha = %.2f, Q_0 = %.2f';
            models(idx).do_include = true;

            idx = idx + 1;
            models(idx).which_structures = 'flat_collins'; 
            models(idx).name = 'RL'; % (flat) Q-learning, as in Collins 2016
            models(idx).params_file = fullfile('results', 'fit_params_results_flat_collins_25nstarts_Q0.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\eta = %.2f, \\beta = %.2f, Q_0 = %.2f';
            models(idx).do_include = true;

            idx = idx + 1;
            models(idx).which_structures = [1 1 1 0 0]; 
            models(idx).name = 'M1, M2, M3';
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M3_25nstarts.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 1 0 1 0]; 
            models(idx).name = 'M1, M2, M1''';
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 1 1 0 0]; 
            models(idx).name = 'M1, M2, M3';
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M3_25nstarts_tau.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f, \\tau^2 = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 1 0 1 0]; 
            models(idx).name = 'M1, M2, M1''';
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f, \\tau^2 = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 1 1 0 0]; 
            models(idx).name = 'M1, M2, M3';
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M3_25nstarts_tau_w0.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f, \\tau^2 = %.2f, w_0 = %.2f';
            models(idx).do_include = false;

            idx = idx + 1;
            models(idx).which_structures = [1 1 0 1 0]; 
            models(idx).name = 'M1, M2, M3'; % our causal structure learning model (note M1' is called M3 in the paper)
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f, \\tau^2 = %.2e, w_0 = %.2f';
            models(idx).do_include = true;

            idx = idx + 1;
            models(idx).which_structures = 'Q_learning'; 
            models(idx).name = 'RL + generalization'; % Q learning with generalization, as proposed by reviewer #1
            models(idx).params_file = fullfile('results', 'fit_params_results_q_learning_2_25nstarts.mat');
            models(idx).params_format = '\\eta = %.2f, \\beta = %.2f, Q_0 = %.2f';
            models(idx).params_idx = 1;
            models(idx).do_include = true;

            idx = idx + 1;
            models(idx).which_structures = [1 0 0 0 0]; 
            models(idx).name = 'M1'; % single structure M1
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_separate_25nstarts_tau_w0.mat');
            models(idx).params_idx = 1;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f, \\tau^2 = %.2e, w_0 = %.2f';
            models(idx).do_include = true;

            idx = idx + 1;
            models(idx).which_structures = [0 1 0 0 0]; 
            models(idx).name = 'M2'; % single structure M2
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_separate_25nstarts_tau_w0.mat');
            models(idx).params_idx = 2;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f, \\tau^2 = %.2e, w_0 = %.2f';
            models(idx).do_include = true;

            idx = idx + 1;
            models(idx).which_structures = [0 0 0 1 0]; 
            models(idx).name = 'M3'; % single structure M1' (we call it M3 in paper but M1' here just to confuse ourselves)
            models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_separate_25nstarts_tau_w0.mat');
            models(idx).params_idx = 3;
            models(idx).params_format = '\\sigma^2_w = %.2f, \\beta = %.2f, \\tau^2 = %.2e, w_0 = %.2f';
            models(idx).do_include = true;

            % filter models -- only include some of them
            models = models(logical([models.do_include]));

            [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

            % load & run models, compute log liks, etc
            %
            for i = 1:numel(models)
                [r, p] = get_test_choice_correlations(models(i).params_file, models(i).params_idx, models(i).which_structures);
               
                load(models(i).params_file, 'results', 'results_options');
                params = results(models(i).params_idx).x;
                options = results_options(models(i).params_idx);
                assert(isequal(models(i).which_structures, options.which_structures));

                total_loglik = model_likfun(data, metadata, params, options.which_structures, data.which_rows, false);
                test_loglik = model_likfun(data, metadata, params, options.which_structures, data.which_rows, true);

                bic = results(models(i).params_idx).bic;

                % sumarize for table
                %
                models(i).params = params;
                models(i).params_string = sprintf(models(i).params_format, params);
                models(i).bic = bic;
                models(i).pxp = NaN;
                models(i).total_loglik = total_loglik;
                models(i).test_loglik = test_loglik;
                models(i).r = r;
                models(i).p = p;
                models(i).p_pow10 = ceil(log10(p));
            end
        
            % compute PXP for pilot data
            % first compute BIC for each subject manually (b/c we did fixed effects => have only one set of params & bic's)
            % need this to compute the PXPs
            %
            [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);

            pilot_lmes = []; % log model evidence
            for i = 1:numel(models)
                K = length(params);
                subj_bics = []; % bic for each subject using the shared fixed effects params
                for who = metadata.subjects % only include "good" subjects
                    which_rows = strcmp(data.participant, who);
                    N = sum(which_rows);

                    subj_loglik = model_likfun(data, metadata, models(i).params, models(i).which_structures, which_rows, false); % from mfit_optimize.m
                    subj_bic = K*log(N) - 2*subj_loglik;
                    subj_bics = [subj_bics; subj_bic];
                end
                pilot_lmes = [pilot_lmes, -0.5 * subj_bics];
                fprintf('Model %d bics (pilot)\n', i);
                disp(subj_bics);
            end
            assert(size(pilot_lmes, 1) == numel(metadata.subjects)); % rows = subjects
            assert(size(pilot_lmes, 2) == numel(models)); % cols = models

            [pilot_alpha,pilot_exp_r,pilot_xp,pilot_pxp,pilot_bor] = bms(pilot_lmes);
            disp('Pilot PXP');
            disp(pilot_pxp);


            % compute PXP for fmri data
            % notice that we don't account for # of parameters here (b/c we fit using pilot data)
            %
            [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

            lmes = []; % log model evidence
            for i = 1:numel(models)
                K = length(params);
                subj_lmes = [];
                for who = metadata.subjects % only include "good" subjects
                    which_rows = strcmp(data.participant, who);
                    N = sum(which_rows);

                    subj_loglik = model_likfun(data, metadata, models(i).params, models(i).which_structures, which_rows, false); % from mfit_optimize.m
                    subj_lmes = [subj_lmes; subj_loglik];
                end
                lmes = [lmes, subj_lmes];
                fprintf('Model %d lmes (fmri)\n', i);
                disp(subj_lmes);
            end
            assert(size(lmes, 1) == numel(metadata.subjects)); % rows = subjects
            assert(size(lmes, 2) == numel(models)); % cols = models

            [alpha,exp_r,xp,pxp,bor] = bms(lmes);
            disp('fMRI PXP');
            disp(pxp);
    
            save(cached_file);
        end

        % Output table
        %
        disp('Model & params & BIC & pPXP & fPXP & Pearson''s r\\');
        for i = 1:numel(models)
            models(i).pilot_pxp = pilot_pxp(i);
            models(i).pxp = pxp(i);
            if models(i).p > 0.0001
                p_string = sprintf('p = %.4f', models(i).p);
            else
                p_string = sprintf('p < 10^{%.0f}', ceil(log10(models(i).p)));
            end
            fprintf('%s & $%s$ & %.0f & %.4f & %.4f & $r = %.2f, %s$ \\\\ \n', ...
                models(i).name, ...
                models(i).params_string, ...
                models(i).bic, ...
                models(i).pilot_pxp, ...
                models(i).pxp, ...
                models(i).r, ...
                p_string);
        end

    case 'KL_structures_154_vs_156'
        
        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results.mat'), 1, [1 1 1 0 ]);
        which_rows = data.which_rows & data.isTrain;
        KL_154 = simulated.surprise(which_rows);

        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_reviewer2.mat'), 1, [1 1 0 1 0]);
        which_rows = data.which_rows & data.isTrain;
        KL_156 = simulated.surprise(which_rows);

        figure;
        %scatter(KL_154, KL_156);
        %lsline;
        plotregression(KL_154, KL_156);
        xlabel('KL structures M1,M2,M3');
        ylabel('KL structures M1,M2,M1''');

        [r, p] = corrcoef(KL_154, KL_156)

    case 'KL_162_vs_163'
        
        [data, metadata, simulated] = simulate_subjects_helper(true, 'results/fit_params_results_simple_collins_25nstarts_0-10alpha_Q0.mat', 1, 'simple_collins');
        which_rows = data.which_rows & data.isTrain;
        KL_162 = simulated.surprise_Zc_given_c(which_rows) + simulated.surprise_Zs_given_s(which_rows);

        [data, metadata, simulated] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, [1 1 0 1 0]);
        which_rows = data.which_rows & data.isTrain;
        KL_163 = simulated.surprise(which_rows);

        figure;
        plotregression(KL_162, KL_163);
        xlabel('KL clusters (Collins)');
        ylabel('KL structures (M1,M2,M1'')');

        [r, p] = corrcoef(KL_162, KL_163)


    case 'tab:glms' % ccnl_bic_bms
        %glms = [161:165 169];
        % 161 = Collins 2016 BUT w/o orth (SPE vs. FPE)
        % 162 = Collins like ours (KL_clusters and PEs), no orth
        % 163 = ours w/ summed KL weights (KL_structures & KL_weights)
        % 164 = ours w/ summed KL weights, weighted by posterior
        % 165 = ours w/ MAP KL_weights
        % 169 = Collins 2016 as is (w/ orth) (SPE vs. FPE)
        % 177 = ours, value & PE (no orth), suggested by reviewer
        %

        %glms = [161 162 163 164 165 177];
        %glms = [163 164 165 177 162 161]; % ordererd as in the paper

        load_cached_values = false;
        cached_file = fullfile('results', 'show_figure_tab_glms.mat');

        if load_cached_values
            load(cached_file);
        else

            idx = 0;

            idx = idx + 1;
            glm(idx).glmodel = 163;
            glm(idx).name = 'GLM 1';
            glm(idx).model = 'M1, M2, M3';
            glm(idx).pmods = '$KL_{structures}$, $KL_{weights}$ (sum)';

            idx = idx + 1;
            glm(idx).glmodel = 164;
            glm(idx).name = 'GLM 2';
            glm(idx).model = 'M1, M2, M3';
            glm(idx).pmods = '$KL_{structures}$, $KL_{weights}$ (weighted sum)';

            idx = idx + 1;
            glm(idx).glmodel = 165;
            glm(idx).name = 'GLM 3';
            glm(idx).model = 'M1, M2, M3';
            glm(idx).pmods = '$KL_{structures}$, $KL_{weights}$ (MAP)';

            idx = idx + 1;
            glm(idx).glmodel = 177;
            glm(idx).name = 'GLM 4';
            glm(idx).model = 'M1, M2, M3';
            glm(idx).pmods = '$V_n$, $PE = r_n - V_n$';

            idx = idx + 1;
            glm(idx).glmodel = 162;
            glm(idx).name = 'GLM 5';
            glm(idx).model = 'RL + clustering';
            glm(idx).pmods = '$KL_{clusters}$, $PE$';

            idx = idx + 1;
            glm(idx).glmodel = 161;
            glm(idx).name = 'GLM 6';
            glm(idx).model = 'RL + clustering$, $RL';
            glm(idx).pmods = 'SPE, FPE';


            bics = [];
            for i = 1:numel(glm)
                bic = ccnl_bic(context_expt(), glm(i).glmodel, 'masks/mask.nii', getGoodSubjects());
                glm(i).bic = bic;
                bics = [bics bic];
            end

            lme = -0.5 * bics;

            [alpha, exp_r, xp, pxp, bor] = bms(lme);

            save(cached_file);
        end

        % output table
        %
        disp('GLM & model & pmods & PXP\\');
        for i = 1:numel(glm)
            glm(i).pxp = pxp(i);
            fprintf('%s & %s & %s & %.4f \\\\ \n', ...
                glm(i).name, ...
                glm(i).model, ...
                glm(i).pmods, ...
                glm(i).pxp);
        end

        pxp

       

    case 'fig:behavior'
        
        %
        % Figure 3A: Posterior probabilities over structures in each condition
        %
        
        figure;
        %set(handle, 'Position', [500, 500, 450, 200])
     
        % M1, M2, M1'
        which_structures = logical([1 1 0 1 0]);
        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat'), 1, [1 1 0 1 0]);
        %[data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_reviewer2.mat'), 1, which_structures);
        
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
        legend({'M1', 'M2', 'M1'''}, 'Position', [0.15 0.3 1 1]);
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        
        text(0.1, 1.25, 'A', 'FontSize', 20, 'FontWeight', 'bold')

        
        %
        % Figure 3B: Choice probabilities on test trials for model vs. humans
        %
        
        subplot(2, 1, 2);

        plot_behavior_helper(data, metadata, simulated);


    case 'RTs_mod'
        % show that modulatory is slower than irr and add, controlling for # of Collins clusters => attention matters => M1,M2,M1'
        % doesn't apply -- hard to control for # of clusters, e.g. if x1c1 came last in the test trials, you already have 3 clusters for each dimension
        %
        [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

        figure;

        subplot(2,1,1);

        which_irr = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'irrelevant') & data.cueId == 2 & data.contextId == 2;
        which_add = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'additive') & data.cueId == 2 & data.contextId == 2;
        which_mod = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'modulatory') & data.cueId == 0 & data.contextId == 0;

        rt_irr_add = data.response.rt(which_irr | which_add);
        rt_irr = data.response.rt(which_irr);
        rt_mod = data.response.rt(which_mod);

        RT_means = [mean(rt_irr) mean(rt_mod)];
        RT_sems = [sem(rt_irr) sem(rt_mod)];

        h = bar(RT_means);
        hold on;
        xs = h(1).XData;
        errorbar(xs, RT_means, RT_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        hold off;
        xticklabels({'irr x3c3', 'mod x1c1'});
        ylim([1 2]);
        ylabel('RT (s)')

        [h, p, ci, stats] = ttest2(rt_irr_add, rt_mod)

        % across subjects -- compare within-subject rt_mod - mean(rt_add, rt_irr) with 0
        %
        subplot(2,1,2);

        mod = [];
        irr = [];
        for who = metadata.subjects
            which_irr = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'irrelevant') & data.cueId == 2 & data.contextId == 2 & strcmp(data.participant, who);
            which_add = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'additive') & data.cueId == 2 & data.contextId == 2 & strcmp(data.participant, who);
            which_mod = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'modulatory') & data.cueId == 0 & data.contextId == 0 & strcmp(data.participant, who);

            rt_irr_add = data.response.rt(which_irr | which_add);
            rt_irr = data.response.rt(which_irr);
            rt_mod = data.response.rt(which_mod);

            mod = [mod mean(rt_mod)];
            irr = [irr mean(rt_irr)];
        end

        RT_means = [mean(mod - irr)];
        RT_sems = [sem(mod - irr)];

        h = bar(RT_means);
        hold on;
        xs = h(1).XData;
        errorbar(xs, RT_means, RT_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        hold off;
        xticklabels({'irr/add x3c3', 'mod x1c1'});
        ylabel('RT (s)')
        title('done properly');

        [h, p, ci, stats] = ttest(mod - irr)

        save shit.mat;
         

    case 'RTs_add'
        % show that additive == irrelevant in terms of feature attention; it's just that 
        % processing text takes longer

        [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

        figure;

        subplot(2,1,1);

        which_irr_old = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'irrelevant') & data.cueId == 0 & data.contextId == 2;
        which_irr_new = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'irrelevant') & data.cueId == 2 & data.contextId == 0;

        which_add_old = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'additive') & data.cueId == 2 & data.contextId == 0;
        which_add_new = data.which_rows & ~data.isTrain & strcmp(data.contextRole, 'additive') & data.cueId == 0 & data.contextId == 2;

        rt_irr_old = data.response.rt(which_irr_old);
        rt_irr_new = data.response.rt(which_irr_new);
        rt_add_old = data.response.rt(which_add_old);
        rt_add_new = data.response.rt(which_add_new);

        RT_means = [mean(rt_irr_old) mean(rt_irr_new) mean(rt_add_old) mean(rt_add_new)];
        RT_sems = [sem(rt_irr_old) sem(rt_irr_new) sem(rt_add_old) sem(rt_add_new)];

        h = bar(RT_means);
        hold on;
        xs = h(1).XData;
        errorbar(xs, RT_means, RT_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        hold off;
        xticklabels({'irr old', 'irr new', 'add old', 'add new'});

        save shit.mat;

        %[h, p, ci, stats] = ttest2(rt_irr, rt_add)


        subplot(2,1,2);

        rt_irr_diff = rt_irr_new - rt_irr_old;
        rt_add_diff = rt_add_new - rt_add_old;

        RT_means = [mean(rt_irr_diff) mean(rt_irr_diff)];
        RT_sems = [sem(rt_irr_diff) sem(rt_irr_diff)];

        h = bar(RT_means);
        hold on;
        xs = h(1).XData;
        errorbar(xs, RT_means, RT_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        hold off;
        xticklabels({'irr diff', 'add diff'});


    case 'RTs_condition'
        % RTs across conditions

        figure;
        %set(handle, 'Position', [500, 500, 450, 200])
     
        % M1, M2, M1'
        which_structures = logical([1 1 0 1 0]);
        %[data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat'), 1, which_structures);
        [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

        %
        % training trials
        %

        subplot(2,1,1);

        RT_means = [];
        RT_sems = [];
        % TODO random effects? LME? repeated measures?
        rt = [];
        g = [];
        i = 0;
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.isTrain & strcmp(data.contextRole, condition) & ~data.timeout;
          
            RTs = data.response.rt(which_rows);
            RT_means = [RT_means, mean(RTs)];
            RT_sems = [RT_sems, std(RTs)/sqrt(length(RTs))];

            i = i + 1;
            g = [g; repmat(i, numel(RTs), 1)];
            rt = [rt; RTs];
        end

        save shit.mat;

        h = bar(RT_means);
        hold on;
        xs = h(1).XData;
        errorbar(xs, RT_means, RT_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        hold off;
        xticklabels({'irr', 'mod', 'add'});
        title('training');
        ylim([1 1.5]);
        ylabel('RT (s)');

        %[p, anovatab, stats] = anova1(rt, g)

        mean(rt(g == 1))
        mean(rt(g == 3))
        [h, p, ci, stats] = ttest2(rt(g == 1), rt(g == 3))


        %
        % test trials
        %

        subplot(2,1,2);

        RT_means = [];
        RT_sems = [];
        % TODO random effects? LME? repeated measures?
        rt = [];
        g = [];
        i = 0;
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.isTest & strcmp(data.contextRole, condition) & ~data.timeout;
          
            RTs = data.response.rt(which_rows);
            RT_means = [RT_means, mean(RTs)];
            RT_sems = [RT_sems, std(RTs)/sqrt(length(RTs))];

            i = i + 1;
            g = [g; repmat(i, numel(RTs), 1)];
            rt = [rt; RTs];
        end

        h = bar(RT_means);
        hold on;
        xs = h(1).XData;
        errorbar(xs, RT_means, RT_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
        hold off;
        xticklabels({'irr', 'mod', 'add'});
        ylim([1 2.2]);
        title('test');
        ylabel('RT (s)');

        %[p, anovatab, stats] = anova1(rt, g)

        mean(rt(g == 1))
        mean(rt(g == 3))
        [h, p, ci, stats] = ttest2(rt(g == 1), rt(g == 3))

    case 'collins_curves'
        %
        % Learning curves, model vs. subjects
        % 
        figure('pos', [100 100 693+20 320] * 3/4);

        %
        % Collins model
        %
     
        % pilot 
        [data, metadata, simulated] = simulate_subjects_helper(false, fullfile('results', 'fit_params_results_simple_collins_25nstarts_0-10alpha_Q0.mat'), 1, 'simple_collins');
        subplot(2,2,1);
        plot_curves_helper(data, metadata, simulated);
        title('Collins: Behavioral pilot');
        text(-4, 1.05, 'A', 'FontSize', 20, 'FontWeight', 'bold')

        % fmri
        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_simple_collins_25nstarts_0-10alpha_Q0.mat'), 1, 'simple_collins');
        subplot(2,2,2);
        plot_curves_helper(data, metadata, simulated);
        title('Collins: fMRI');
        text(-4, 1.05, 'B', 'FontSize', 20, 'FontWeight', 'bold')

        %
        % Flat model
        %
     
        % pilot 
        [data, metadata, simulated] = simulate_subjects_helper(false, fullfile('results', 'fit_params_results_flat_collins_5nstarts.mat'), 1, 'flat_collins');
        subplot(2,2,3);
        plot_curves_helper(data, metadata, simulated);
        title('Flat: Behavioral pilot');
        text(-4, 1.05, 'A', 'FontSize', 20, 'FontWeight', 'bold')

        % fmri
        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_flat_collins_5nstarts.mat'), 1, 'flat_collins');
        subplot(2,2,4);
        plot_curves_helper(data, metadata, simulated);
        title('Flat: fMRI');
        text(-4, 1.05, 'B', 'FontSize', 20, 'FontWeight', 'bold')


    case 'collins_behavior'

        figure
        
        % Collins & frank model
        %

        subplot(2,1,1);

        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_simple_collins_25nstarts_0-10alpha_Q0.mat'), 1, 'simple_collins');

        plot_behavior_helper(data, metadata, simulated);
        title('Collins');

        % Flat Q learning
        %

        subplot(2,1,2);

        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_flat_collins_5nstarts.mat'), 1, 'flat_collins');

        plot_behavior_helper(data, metadata, simulated);
        title('Flat RL');

    case 'collins_posteriors'
        
        
        figure;
        %set(handle, 'Position', [500, 500, 450, 200])
     
        % M1, M2, M1'
        %which_structures = logical([1 1 0 1 0]);
        %[data, metadata, simulated] = simulate_subjects_helper();        
        which_structures = 'simple_collins';
        [data, metadata, simulated] = simulate_subjects_helper(true, fullfile('results', 'fit_params_results_simple_collins_5nstarts.mat'), 1, which_structures);
        

        i = 0;
        P_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.newTrialId == 20 & strcmp(data.contextRole, condition);
            
            P = simulated.posteriors_Zc_given_C(:,:,which_rows);
            P = mean(P, 3);
            P = mean(P, 2)';
            P = P(1:3);
            P_means = [P_means; P];

            i = i + 1;
        end
        subplot(4,1,1);
        bar(P_means);
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('P(Z_c) ');
        legend({'Z_c1', 'Z_c2', 'Z_c3'});
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        title('context cluster popularities after training');
        

        
        i = 0;
        P_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.newTrialId == 24 & strcmp(data.contextRole, condition);
            
            P = simulated.posteriors_Zc_given_C(:,:,which_rows);
            P = mean(P, 3);
            P = mean(P, 2)';
            P = P(1:3);
            P_means = [P_means; P];

            i = i + 1;
        end
        subplot(4,1,2);
        bar(P_means);
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('P(Z_c) ');
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        title('context cluster popularities after test');




        i = 0;
        P_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.newTrialId == 20 & strcmp(data.contextRole, condition);
            
            P = simulated.posteriors_Zs_given_S(:,:,which_rows);
            P = mean(P, 3);
            P = mean(P, 2)';
            P = P(1:3);
            P_means = [P_means; P];

            i = i + 1;
        end
        subplot(4,1,3);
        bar(P_means);
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        legend({'Z_s1', 'Z_s2', 'Z_s3'});
        ylabel('P(Z_s) ');
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        title('cue cluster popularities after training');
        

        
        i = 0;
        P_means = [];
        for condition = metadata.contextRoles
            which_rows = data.which_rows & data.newTrialId == 24 & strcmp(data.contextRole, condition);
            
            P = simulated.posteriors_Zs_given_S(:,:,which_rows);
            P = mean(P, 3);
            P = mean(P, 2)';
            P = P(1:3);
            P_means = [P_means; P];

            i = i + 1;
        end
        subplot(4,1,4);
        bar(P_means);
        xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
        ylabel('P(Z_s) ');
        ylim([0 1.1]);
        set(gca,'fontsize',13);
        title('cue cluster popularities after test');




    case '_q_behavior'

        % behavioral plot for simple Q learning as suggested by reviewer 1
        %

        models(1).which_structures = 'simple_Q'; 
        models(1).name = 'Q learning';
        models(1).params_file = fullfile('results', 'fit_params_results_simple_q.mat');
        models(1).params_format = '\\beta = %.4f';
        models(1).params_idx = 1;

        models(2).which_structures = 'Q_learning'; 
        models(2).name = 'Q learning 2';
        models(2).params_file = fullfile('results', 'fit_params_results_q_learning.mat');
        models(2).params_format = '\\alpha = %.4f, \\beta = %.4f';
        models(2).params_idx = 1;

        for i=1:numel(models)

            % plot learning curves
            %
            figure;

            % pilot 
            [data, metadata, simulated] = simulate_subjects_helper(false, models(i).params_file, models(i).params_idx, models(i).which_structures);
            subplot(2,2,1);
            plot_curves_helper(data, metadata, simulated);
            title('Behavioral pilot');
            text(-4, 1.05, 'A', 'FontSize', 20, 'FontWeight', 'bold')
               
            % fmri
            [data, metadata, simulated] = simulate_subjects_helper(true, models(i).params_file, models(i).params_idx, models(i).which_structures);
            subplot(2,2,2);
            plot_curves_helper(data, metadata, simulated);
            title('fMRI');
            text(-4, 1.05, 'B', 'FontSize', 20, 'FontWeight', 'bold')

            % plot test choices
            %
            [data, metadata, simulated] = simulate_subjects_helper(true, models(i).params_file, models(i).params_idx, models(i).which_structures);
            subplot(2,1,2);
            
            plot_behavior_helper(data, metadata, simulated);

            title(models(i).name);
        end
        
    case 'searchlight_posterior'
        bspmview('rdms/betas_smooth/searchlight_tmap_posterior_feedback_onset.nii', 'masks/mean.nii');
        
    case 'searchlight_prior'
        bspmview('rdms/betas_smooth/searchlight_tmap_prior_trial_onset.nii', 'masks/mean.nii');
        
    case 'KL_structures'
        ccnl_view(context_expt(), 154, 'KL_structures');
        
    case 'KL_weights'
        ccnl_view(context_expt(), 154, 'KL_weights');
        
    case 'KL_structures - KL_weights'
        ccnl_view(context_expt(), 154, 'KL_structures - KL_weights');

    case 'KL_weights - KL_structures'
        ccnl_view(context_expt(), 154, 'KL_weights - KL_structures');
        
    case 'fig:glm163'
        figure('pos', [100 100 653 492]);

        p = panel();
        p.pack(1, 3);
        p.de.margin = 0;

        fontsize = 12;
        
        p(1, 1).select();
        imshow('images/glm163_KL_structures.png'); 
        title('KL_{structures}', 'FontSize', fontsize);

        p(1, 2).select();
        imshow('images/glm163_KL_weights.png'); 
        title('KL_{weights}', 'FontSize', fontsize);

        p(1, 3).select();
        imshow('images/glm163_KL_structures-KL_weights.png'); 
        title('KL_{structures} - KL_{weights}', 'FontSize', fontsize);

        p.export('../JNeuro manuscript/figures/glm163.pdf');


    case 'fig:fmri-results'
    %case 'glm154'
        figure('pos', [100 100 693+20 492]);
        
        fontsize = 12;
        
        %
        % Top left: structures 
        %

        subplot(2, 2, 1);
        
        imshow('images-old/KL_structures_pos.png'); % from GLM 154
        title('Causal structure update', 'FontSize', fontsize);
       
        %
        % Top right: weights 
        %
        
        subplot(2, 2, 2);
        
        imshow('images-old/KL_weights_pos.png'); % from GLM 154
        title('Associative weights update', 'FontSize', fontsize);
        
        %
        % Bottom left: contrast
        %
        
        subplot(2, 2, 3);
        
        imshow('images-old/KL_structures-KL_weights.png'); % from GLM 154
        title({'Causal structure update >'; 'associative weights update'}, 'FontSize', fontsize);
        
        %
        % Bottom right: KL structures ~ test choice log likelihood
        %

        subplot(2, 2, 4);
        
        load results/betas_to_behavior_glm154_KL_structures_sphere_KL_structures.mat
        which = 1:numel(region);
        %assert(numel(region) == 14);
        %which = [1:7 12 14]; % exclude motor and visual areas
        warning('this is shit -- result did not reproduce');
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
        

    case 'fig:rsa'
        
        figure('pos', [100 100 693+20 492]);
        
        fontsize = 12;
        
        %
        % Top left: structures 
        %

        subplot(1, 2, 1);
        
        imshow('images-old/searchlight_prior.png');
        title('Causal structure prior', 'FontSize', fontsize);
        
        g = subplot(1, 2, 2);
        
        imshow('images-old/searchlight_posterior.png');
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

        
        
        







        
    case '_Figure_4A_Neuron_DEPRECATED'
        ccnl_view(context_expt(), 123, 'surprise');
        
    case '_Figure_4C_Neuron_DEPRECATED'
        ccnl_view(context_expt(), 148, 'KL_weights');
        

    case '_Figure_4A_stats_Neuron_DEPRECATED'        
        
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
                
        
        
        
        
    case '_Figure_4_Neuron_DEPRECATED'
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
        
        
        
        
        
        
        
    case '_Figure_4B_DEPRECATED'     
        
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


    case '_Figure_3_stats_DEPRECATED'        
        
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
       
        % TODO
        assert(false, 'THESE ARE DEPRECATED -- use the stats from plot_behavior.m; you have to change which_structures e.g. to [1 0 0 0] for M1 only');
        
    end 


    save(fullfile('results', 'show_figure_last_run.mat'));

end % show_figure()









% helper f'n to plot behavioral results
%
function plot_behavior_helper(data, metadata, simulated)

    % Choice probabilities for model
    %
    model_means = [];
    model_sems = [];
    for condition = metadata.contextRoles
        which_rows = data.which_rows & ~data.isTrain & strcmp(data.contextRole, condition);
        
        x1c1 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 0);
        x1c3 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 2);
        x3c1 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 0);
        x3c3 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 2);

        model_means = [model_means; mean(x1c1) mean(x1c3) mean(x3c1) mean(x3c3)];
        model_sems = [model_sems; std(x1c1)/sqrt(numel(x1c1)) std(x1c3)/sqrt(numel(x1c3)) std(x3c1)/sqrt(numel(x3c1)) std(x3c3)/sqrt(numel(x3c3))];
    end

    % Plot model choices probabilities
    %
    h = bar(model_means);

    hold on;

    if nargin == 2 & ~isempty(model_sems)
        xs = sort([h(1).XData + h(1).XOffset, ...
                   h(2).XData + h(2).XOffset, ...
                   h(3).XData + h(3).XOffset, ...
                   h(4).XData + h(4).XOffset]);
        model_means = model_means'; model_means = model_means(:);
        model_sems = model_sems'; model_sems = model_sems(:);        
        errorbar(xs, model_means, model_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    end

    % plot both pilot and fmri subjects' choices
    %
    group_x_offs = [-0.02, 0.02];
    group_color = {[0.5 0.5 0.5], [0 0 0]};
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
end



% helper for plotting learning curves
%
function plot_curves_helper(data, metadata, simulated, which_rows)
    sem = @(x) std(x) / sqrt(length(x));

    if nargin < 4 || isempty(which_rows)
        which_rows = data.which_rows;
    end

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
            which = which_rows & data.isTrain & data.trialId == n & strcmp(data.participant, who);
            %assert(sum(which) == 9);
            
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

    plot(model_correct, '.-', 'LineWidth', 2); % == mean(human_correct_all_runs)
    %errorbar(1:metadata.trainingTrialsPerRun, model_correct, model_correct_sem, 'o-', 'LineWidth', 2);
    hold on;
    plot(human_correct, '.-', 'LineWidth', 2); % == mean(model_correct_all_runs)
    %errorbar(1:metadata.trainingTrialsPerRun, human_correct, human_correct_sem, 'o-', 'LineWidth', 2);
    hold off;
    legend({'model', 'subjects'}, 'Position', [0.30 0.05 1 1]);
    %title('Average per-trial accuracy');
    xlabel('trial #');
    ylabel('accuracy');
    set(gca,'fontsize',13);
end




% correlate model choices with subject choices on the test trials
%
function [r, p] = get_test_choice_correlations(params_file, params_idx, which_structures)
    utils; % include some nifty lambdas

    [data, metadata, simulated] = simulate_subjects_helper(true, params_file, params_idx, which_structures);

    %
    % Choice probabilities in test phase for SUBJECTS
    %

    Ms = [];
    SEMs = [];
    for context = metadata.contextRoles
        which = data.which_rows & data.isTrain == 0 & strcmp(data.contextRole, context);

        x1c1 = strcmp(data.response.keys(which & data.cueId == 0 & data.contextId == 0), 'left');
        x1c2 = strcmp(data.response.keys(which & data.cueId == 0 & data.contextId == 2), 'left');
        x2c1 = strcmp(data.response.keys(which & data.cueId == 2 & data.contextId == 0), 'left');
        x2c2 = strcmp(data.response.keys(which & data.cueId == 2 & data.contextId == 2), 'left');

    %    M = mean([x1c1 x1c2 x2c1 x2c2]);
    %    SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
        M = get_means(x1c1, x1c2, x2c1, x2c2);
        SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
        Ms = [Ms; M];
        SEMs = [SEMs; SEM];
    end

    subject_Ms = Ms; % for stats

    %
    % TRUE Choice probabilities in test phase for MODEL
    %

    Ms = [];
    SEMs = [];
    for context = metadata.contextRoles
        which = data.which_rows & data.isTrain == 0 & strcmp(data.contextRole, context);

        x1c1 = simulated.pred(which & data.cueId == 0 & data.contextId == 0);
        x1c2 = simulated.pred(which & data.cueId == 0 & data.contextId == 2);
        x2c1 = simulated.pred(which & data.cueId == 2 & data.contextId == 0);
        x2c2 = simulated.pred(which & data.cueId == 2 & data.contextId == 2);

        %M = mean([x1c1 x1c2 x2c1 x2c2]);
        %SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
        M = get_means(x1c1, x1c2, x2c1, x2c2);
        SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
        Ms = [Ms; M];
        SEMs = [SEMs; SEM];
    end

    model_Ms = Ms; % for stats

    % correlate average subject choices with model choices 
    %
    [r, p] = corrcoef(subject_Ms(:), model_Ms(:));
    r = r(1,2);
    p = p(1,2);
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
