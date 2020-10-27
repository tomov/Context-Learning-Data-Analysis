
clear all;

[params, which_structures] = model_params('results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat')

[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

%{
load rdms/rdms_searchlight_isl.mat

for i = 1:length(Model)
    if startsWith(Model(i).name, 'MCMC')
        simulated = Model(i).simulated;
        model_name = Model(i).name;
        filename = sprintf('%s_np=1000.mat', model_name);
        save(fullfile('mat', filename), 'simulated', 'model_name', 'params');
    end
end

model_name = 'MCMC_neurath4';    
simulated = simulate_subjects(data, metadata, params, model_name);
filename = sprintf('%s_np=1000.mat', model_name);
save(fullfile('mat', filename), 'simulated', 'model_name', 'params');
%}



model_name = 'MCMC_ideal_w=1';    
simulated = simulate_subjects(data, metadata, params, model_name);
filename = sprintf('%s_np=1000.mat', model_name);
save(fullfile('mat', filename), 'simulated', 'model_name', 'params');

model_name = 'MCMC_reset_w=1';    
simulated = simulate_subjects(data, metadata, params, model_name);
filename = sprintf('%s_np=1000.mat', model_name);
save(fullfile('mat', filename), 'simulated', 'model_name', 'params');

model_name = 'MCMC_neurath_w=1';    
simulated = simulate_subjects(data, metadata, params, model_name);
filename = sprintf('%s_np=1000.mat', model_name);
save(fullfile('mat', filename), 'simulated', 'model_name', 'params');
