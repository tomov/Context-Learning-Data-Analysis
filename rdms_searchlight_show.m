% Create a t-map for a given event and model, based on the pre-computed searchlight results by rdms_searchlight.m in the rdms/ directory
%

clear all;
close all;

% QUESTIONS:
% 1. P(weights) -- comparing trials across runs? see rdms_get_model.m
% NO -- too similar; 2. two-sample t-test of prior vs. posterior -- which one is more similar to the
% given voxel RDM? (we have spearman's rho for all participants)
%

event = 'feedback_onset';% <---------------------------- name of the game
model = 'posterior';% <---------------------------- name of the game
%event = 'trial_onset';
%model = 'prior';
%model = 'ww_prior';
%event = 'feedback_onset';
%model = 'ww_posterior';
%dirname = 'rdms/betas_smooth';
%dirname = 'rdms/M1M2M1_25nstarts_tau_w0';
%event = 'feedback_onset';
%model = 'Q_posteosterior';
%model = 'cond_posterior';

dirname = 'rdms/M1M2M1_4mm'; % <---------------------------- name of the game
%dirname = 'rdms/M1M2M1_4mm';
%dirname = 'rdms/M1M2M1_4mm_nosmooth';



%event = 'feedback_onset';
%model = 'joint_posterior'; % nothing
%model = 'cond_posterior'; % nothing
%model = 'Q_posteosterior'; % IFG
%event = 'trial_onset';
%model = 'joint_prior'; % nothing
%model = 'cond_prior'; % nothing
%model = 'Q_prior'; % nothing 

%dirname = 'rdms/Collins_4mm';

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
Model = rdms_get_model_3(data, metadata, which_rows); % WARNING == TIGHT COUPLING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! sanity check if params are same as in the searchlight file
%Model = rdms_get_model_collins(data, metadata, which_rows); % WARNING == TIGHT COUPLING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! sanity check if params are same as in the searchlight file
model_names = {Model.name};
assert(ismember(model, model_names));

%% initialize an empty t-map
%
[~, V, tmap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
tmap(:) = NaN; % clear
pmap = tmap;
V.fname = fullfile(dirname, ['searchlight_tmap_', model, '_', event, '.nii']); % change immediately!

maxt = 0;

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
            % a two-tailed t-test), so we have to use rho to get it.
            % as a reminder, this was a 1-sample t-test that rho != 0
            %
            df = metadata.N - 1;
            table_T = icdf('T', 1 - table_P / 2, df);
            sign = (table_Rho > 0) + (-1) * (table_Rho < 0);
            table_T = table_T .* sign;
            disp(['Computing t-values manually for ', file]);
        end

		% this was used to get the coordinates of the nosmooth searchlights when we didn't record them
        %
        %if exist(file, 'file')
        %    fprintf('Skipping %s -- already exists\n', file);
        %    continue;
        %end
        %Searchlight = rdms_get_searchlight(data, metadata, data.which_rows & data.isTrain, x, y, z, r, true, use_tmaps, use_nosmooth); % use pregen'd betas, use tmaps, use nosmooth
        %x = [Searchlight.x];
        %y = [Searchlight.y];
        %z = [Searchlight.z];
        %events = {Searchlight.event};
        %save(file, 'table_Rho', 'table_T', 'table_P', 'all_subject_rhos', 'x', 'y', 'z', 'events', 'r', 'idx', 'params', 'which_structures', 'use_tmaps', 'use_nosmooth', 'which_rows');
       
        if exist('events', 'var') % new way
            % when using _nosmooth, sometimes not all searchlights are good
            % => have to keep 
            %
            for j = 1:length(x)
                if ~isequal(event, events{j}), continue; end
                row_idx = j;
                col_idx = find(ismember(model_names, model));
                if col_idx <= size(table_T, 2) % in case we're using a model that hasn't been computed in all batches, e.g. weightsPosterior
                    tmap(x(j), y(j), z(j)) = table_T(row_idx, col_idx);
                end
            end
            clear events; % for next iteration
        else % old way
            % backwards compatibility -- we have all the searchlights when using smoothing; back then 
            % we didn't 
            %
            for j = 1:length(x)
                row_idx = j + isequal(event, 'feedback_onset') * length(x);
                col_idx = find(ismember(model_names, model));
                if col_idx <= size(table_T, 2) % in case we're using a model that hasn't been computed in all batches, e.g. weightsPosterior
                    tmap(x(j), y(j), z(j)) = table_T(row_idx, col_idx);
                    pmap(x(j), y(j), z(j)) = table_P(row_idx, col_idx);
                end

                % find the max correlation for displaying
                t =   table_T(row_idx, col_idx);
                if t > maxt
                    maxt = t;
                    maxt_idx = j;
                    maxt_file = file;
                    maxt_x = x(j);
                    maxt_y = y(j);
                    maxt_z = z(j);
                end
            end
        end
    end
end

% FDR
%
%{
alpha = 0.05;
p = pmap(~isnan(pmap));
fdr = mafdr(p, 'BHFDR', true);
pmap(~isnan(pmap)) = fdr;
tmap(pmap >= alpha) = NaN;
%}

spm_write_vol(V, tmap);

%% visualize
%
EXPT = context_expt();
struc = fullfile('masks','mean.nii');

disp(V.fname)

bspmview(V.fname, struc);
