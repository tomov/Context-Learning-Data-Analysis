%clear all;

%CV_expts;

r = 10; % mm
masks = isl_create_masks(false, r);

glmodels = 194:197;
goodSubjs = getGoodSubjects();

for m = 1:length(masks)
    maskfile = masks{m};

    [~,mask,~] = get_mask_format_helper(maskfile);
    inds = find(mask);

    for g = 1:length(glmodels)
        glmodel = glmodels(g);

        for s = 1:length(goodSubjs)
            subj = goodSubjs(s);

            for k = 1:3

                fprintf('m=%d g=%d s=%d k=%d\n', m, g, s, k);

                test_modeldir = fullfile(EXPT_test{k}.modeldir, ['model',num2str(glmodel)],['subj',num2str(subj)]);
                train_modeldir = fullfile(EXPT_train{k}.modeldir, ['model',num2str(glmodel)],['subj',num2str(subj)]);

                load(fullfile(test_modeldir, 'SPM.mat'));
                SPM_test = SPM;

                load(fullfile(train_modeldir, 'SPM.mat'));
                SPM_train = SPM;

                Y_test = spm_data_read(SPM_test.xY.VY, inds);
                KWY_test = spm_filter(SPM_test.xX.K, SPM_test.xX.W * Y_test);

                cdir = pwd;
                cd(train_modeldir);
                B_train = spm_data_read(SPM_train.Vbeta, inds);
                cd(test_modeldir);
                B_test = spm_data_read(SPM_test.Vbeta, inds);
                cd(cdir);

                res = KWY_test - SPM_test.xX.xKXs.X * B_test; % for sanity
                res_CV = KWY_test - SPM_test.xX.xKXs.X * B_train;

                res2 = ccnl_get_residuals(EXPT_test{k}, glmodel, maskfile, subj, true, true);
                res2 = res2{1};
                assert(immse(res, res2) < 1e-9);

                ResMS_CV = sum(res_CV.^2, 1) / SPM_test.xX.trRV;
                mses{m}(s,g,k) = mean(ResMS_CV);

                % sanity checks
                %
                %{
                KWX = spm_filter(SPM_test.xX.K, SPM_test.xX.W * SPM_test.xX.X);
                immse(KWX, SPM_test.xX.xKXs.X)

                beta = SPM_test.xX.pKX * KWY_test;
                b1 = ccnl_get_beta(EXPT_test{k}, glmodel, 'M1', maskfile, subj);
                b2 = ccnl_get_beta(EXPT_test{k}, glmodel, 'M1', mask, subj);
                b3 = ccnl_get_beta(EXPT_test{k}, glmodel, 'M1', inds, subj);
                assert(immse(b1, beta(1,:)) < 1e-9);
                assert(immse(b2, beta(1,:)) < 1e-9);
                assert(immse(b3, beta(1,:)) < 1e-9);
                assert(immse(beta, B_test) < 1e-9);

                %res0 = KWY_test - KWX * B_test; % for sanity

                res3 = spm_sp('r',SPM_test.xX.xKXs,KWY_test);
                assert(immse(res3, res2) < 1e-9);
                %}
            end 
        end
    end

    mse{m} = mean(mses{m}, 3);

end

save('mat/CV.mat', '-v7.3');

plot_CV;
