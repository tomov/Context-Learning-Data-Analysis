% Create a t-map for a given event and model, based on the pre-computed searchlight results by rdms_searchlight.m in the rdms/ directory
%

clear all;
close all;

% QUESTIONS:
% 1. P(weights) -- comparing trials across runs? see rdms_get_model.m
% NO -- too similar; 2. two-sample t-test of prior vs. posterior -- which one is more similar to the
% given voxel RDM? (we have spearman's rho for all participants)
%

%event = 'feedback_onset';
%model = 'ww_Sigma_posterior';
event = 'feedback_onset';
%dirname = 'rdms/betas_smooth';
dirname = 'rdms';

% posterior @ feedback_onset -- bilateral AG, bilateral dlPFC, IT, visual... :(
% prior @ feedback_onset -- same
% prior @ trial_onset -- vlPFC
% values @ feedback_onset - bilateral AG, right dlPFC
% value @ feedback_onset - right IT
% newValue @ feedback_onset - right IT, right AG !!!
% PEs - bilateral AG
% 
% weightsAfter @ feedback_onset -- nothing :(
% weightsBefore @ trial_onset -- 

% logPosterior @ feedback_onset -- bilateral AG, dlPFC
% logPrior @ feedback_onset -- bilateral AG, lateral PFC, visual
% logPrior @ trial_onset -- right IT, med PFC / BA 10
% logValuesSquared @ feedback_onset -- right AG

assert(ismember(event, {'trial_onset', 'feedback_onset'}));

%% get model names
%
[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());
which_rows = data.which_rows & data.isTrain; % Look at training trials only
Model = rdms_get_model_3(data, metadata, which_rows);
model_names = {Model.name};
assert(ismember(model, model_names));

%% initialize an empty t-map
%
[~, V, tmap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
tmap(:) = NaN; % clear
V.fname = fullfile(dirname, ['searchlight_tmap_', model, '_', event, '.nii']); % change immediately!

%% load all the searchlight results
%
files = dir(dirname);
for i = 1:length(files)
    file = files(i).name;
    if startsWith(file, 'searchlight_') && endsWith(file, '.mat')
        disp(['Loading ', file]);
        table_T = [];
        load(fullfile(dirname, file));
        if isempty(table_T)
            % forgot to keep track of table_T sometimes...
            % you can get the t-values from the p-values using the fact that:
            % p = (1 - tcdf(t, df)) * 2
            % so inverting:
            % t = icdf('T', 1 - p / 2, df)
            % note that this doesn't tell you anything about the sign (it's
            % a two-tailed t-test), so we have to use rho to get it
            % as a reminder, this was a 1-sample t-test that rho != 0
            %
            df = metadata.N - 1;
            table_T = icdf('T', 1 - table_P / 2, df);
            sign = (table_Rho > 0) + (-1) * (table_Rho < 0);
            table_T = table_T .* sign;
            disp(['Computing t-values manually for ', file]);
        end
        
        for j = 1:length(x)
            row_idx = j + isequal(event, 'feedback_onset') * length(x);
            col_idx = find(ismember(model_names, model));
            if col_idx <= size(table_T, 2) % in case we're using a model that hasn't been computed in all batches, e.g. weightsPosterior
                tmap(x(j), y(j), z(j)) = table_T(row_idx, col_idx);
            end
        end
    end
end

spm_write_vol(V, tmap);

%% visualize
%
EXPT = context_expt();
struc = fullfile(EXPT.modeldir,'mean.nii');

bspmview(V.fname, struc);
