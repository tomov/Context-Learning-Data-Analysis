idx = 0;

idx = idx + 1;
models(idx).which_structures = [1 1 0 1 0]; 
models(idx).name = 'M1, M2, M3'; % our causal structure learning model (note M1' is called M3 in the paper)
models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat');
models(idx).params_idx = 1;
models(idx).params_format = '\\sigma^2_w = %.4f, \\beta = %.4f, \\tau^2 = %.4e, w_0 = %.4f';
models(idx).do_include = false;

idx = idx + 1;
models(idx).which_structures = 'ideal'; 
models(idx).name = 'ideal';
models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat');
models(idx).params_idx = 1;
models(idx).params_format = '\\sigma^2_w = %.4f, \\beta = %.4f, \\tau^2 = %.4e, w_0 = %.4f';
models(idx).do_include = true;


idx = idx + 1;
models(idx).which_structures = 'MCMC_ideal'; 
models(idx).name = 'MCMC_ideal';
models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat');
models(idx).params_idx = 1;
models(idx).params_format = '\\sigma^2_w = %.4f, \\beta = %.4f, \\tau^2 = %.4e, w_0 = %.4f';
models(idx).do_include = true;

idx = idx + 1;
models(idx).which_structures = 'MCMC_reset'; 
models(idx).name = 'MCMC_reset';
models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat');
models(idx).params_idx = 1;
models(idx).params_format = '\\sigma^2_w = %.4f, \\beta = %.4f, \\tau^2 = %.4e, w_0 = %.4f';
models(idx).do_include = true;

idx = idx + 1;
models(idx).which_structures = 'MCMC_neurath'; 
models(idx).name = 'MCMC_neurath';
models(idx).params_file = fullfile('results', 'fit_params_results_M1M2M1_25nstarts_tau_w0.mat');
models(idx).params_idx = 1;
models(idx).params_format = '\\sigma^2_w = %.4f, \\beta = %.4f, \\tau^2 = %.4e, w_0 = %.4f';
models(idx).do_include = true;

% filter models -- only include some of them
models = models(logical([models.do_include]));

% compute PXP for fmri data
% notice that we don't account for # of parameters here (b/c we fit using pilot data)
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

lmes = []; % log model evidence
for i = 1:numel(models)
    load(models(i).params_file);
    params = results.x;
    assert(size(params,1) == 1);
    K = length(params);
    subj_lmes = [];
    for who = metadata.subjects %(1) % only include "good" subjects
        which_rows = strcmp(data.participant, who); % & data.runId == 2;
        N = sum(which_rows);

        subj_loglik = model_likfun(data, metadata, params, models(i).which_structures, which_rows, false); % from mfit_optimize.m
        subj_lmes = [subj_lmes; subj_loglik];
    end
    lmes = [lmes, subj_lmes];
    fprintf('Model %d lmes (fmri)\n', i);
    disp(subj_lmes);
end
%assert(size(lmes, 1) == numel(metadata.subjects)); % rows = subjects
assert(size(lmes, 2) == numel(models)); % cols = models

[alpha,exp_r,xp,pxp,bor] = bms(lmes);
disp('fMRI PXP');
disp(pxp);
