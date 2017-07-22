% Compare RDMs from different ROIs with model RDMs
%

%% Load data and compute first-order RDMs
%

[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());

which_rows = data.which_rows & data.isTrain; % Look at training trials only

% Get the neural RDMs
%
Neural = rdms_get_neural(data, metadata, which_rows, true, true); % use t-maps, use nosmooth
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
    do_LME = false; % we can't do LME here
    rows = [Neural Model];
    cols = rows;
    [table_Rho, table_H, table_T, table_P, all_subject_rhos] = rdms_second_order(metadata, rows, cols, [], do_LME, [], []);
else
    % Neural RDMs vs. Model RDMs
    %
    % It is important that rows = Neural and cols = Model for the visualizations
    % and the linear mixed effects modeling
    %
    do_LME = true;
    %lme_neural_idxs = 1:numel(Neural);
    lme_neural_idxs = [1 2 5 11 12 14 15 18 24 25]; % which ROIs to consider for LME
    lme_model_idxs = [1 3 18 39 41 46]; % which models to consider to LME
    rows = Neural;
    cols = Model;
    [table_Rho, table_H, table_T, table_P, all_subject_rhos, lmes] = rdms_second_order(metadata, rows, cols, control_model_idxs, do_LME, lme_neural_idxs, lme_model_idxs);
end


Rho = array2table(table_Rho, 'RowNames', {rows.name}, 'VariableNames', {cols.name});
H = array2table(table_H, 'RowNames', {rows.name}, 'VariableNames', {cols.name});
P = array2table(table_P, 'RowNames', {rows.name}, 'VariableNames', {cols.name});




%% Visualize full analysis
%
rdms_show_full(rows, cols, table_Rho, table_P, all_vs_all);

%% Show the significant positive correlations
%
%rdms_show_significant(rows, cols, table_Rho, table_P, all_subject_rhos, lmes, do_LME, all_vs_all, true); % with LME
rdms_show_significant(rows, cols, table_Rho, table_P, all_subject_rhos, {}, false, all_vs_all, true); % no LME

%% Chan et al.-like bar plots for their models
%
rdms_show_chan_etal(rows, cols, all_subject_rhos, all_vs_all, control_model_idxs);

%% Interesting visualizations
% show select neural / model pairs that I think are interesting
%
rdms_show_interesting(rows, cols, all_subject_rhos, all_vs_all, control_model_idxs);
