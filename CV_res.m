
% CV with beta series GLM

% see CV.m

clear all;

r = 10; % mm

masks = isl_create_masks(false, r);
masks = masks(1:4);

glmodel = 198;

goodSubjs = getGoodSubjects();

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

[allSubjects, subjdirs, nRuns] = context_getSubjectsDirsAndRuns();

which_structuress = {logical([1 1 0 1 0]), 'MCMC_ideal', 'MCMC_reset', 'MCMC_neurath'};

for m = 1:length(masks)
    maskfile = masks{m};

    for g = 1:length(which_structuress)
        which_structures = which_structuress{g};

        for s = 1:length(goodSubjs)
            subj = goodSubjs(s);

            which_train = data.isTrain & strcmp(data.participant, allSubjects{subj});
            which_test = data.isTest & strcmp(data.participant, allSubjects{subj});

            [~,~,simulated] = simulate_subjects_helper(true, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 1, which_structures, which_train | which_test);

            P = simulated.P(which_train,:);
            X = P(:, [1 2 4]);

            Y = ccnl_get_beta_series(context_expt(), glmodel, subj, 'feedback_onset', maskfile);

            partition_id = ceil([1:size(Y,1)]/60);
            run_id = ceil([1:size(Y,1)]/20);

            Z = nan(size(Y));
            for r = 1:9
                Z(run_id == r, :) = zscore(Y(run_id == r,:), 0, 1);
            end

            for k = 1:3
                fprintf('m=%d g=%d s=%d k=%d\n', m, g, s, k);

                test = partition_id == k;
                train = partition_id ~= k;

                X_test = X(test,:);
                Z_test = Z(test,:);
                X_train = X(train,:);
                Z_train = Y(train,:);

                b = (X_train' * X_train)^(-1) * X_train' * Z_train; % OLS; see test_gp.m in VGDL matlab repo

                Z_pred = X_test * b;

                mses{m}(s,g,k) = immse(Z_pred, Z_test);
                r = zeros(1, size(Z_pred,2));
                for i = 1:size(Z_pred, 2)
                    r(i) = corr(Z_pred(:,i), Z_test(:,i));
                end
                avg_rs{m}(s,g,k) = mean(r);
                r_avgs{m}(s,g,k) = corr(mean(Z_pred,2), mean(Z_test,2));

            end

        end
    end

    mse{m} = mean(mses{m}, 3);
    avg_r{m} = mean(avg_rs{m}, 3);
    r_avg{m} = mean(r_avgs{m}, 3);
end

save('mat/CV_res_zscore.mat', '-v7.3');

CV_plot;
