% run after isl_bms.m

figure;
for i = 1:numel(models)

    subplot(length(models), 1, i);
    clear simulated;
    clear data; 
    clear metadata;
    [data, metadata, simulated] = simulate_subjects_helper(true, models(i).params_file, 1, models(i).which_structures);
    plot_behavior_helper(data, metadata, simulated);

    title(models(i).name, 'interpreter', 'none');
end
