%[params, which_structures] = model_default_params();
params = [0.1249 2.0064];
which_structures = [1 1 1 0];

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

simulated_old = simulate_subjects_old(data, metadata, params, which_structures);

simulated = simulate_subjects(data, metadata, params, which_structures);

assert(abs(sum(sum(simulated_old.valuess - simulated.valuess(:,1:4)))) < 1e-6);
assert(abs(sum(sum(simulated_old.P - simulated.P(:,1:4)))) < 1e-6);
assert(abs(sum(sum(simulated_old.pred - simulated.pred))) < 1e-6);
assert(abs(sum(sum(simulated_old.values - simulated.values))) < 1e-6);
