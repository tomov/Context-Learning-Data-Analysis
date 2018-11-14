% Create a t-map for a given event and model, based on the pre-computed searchlight results by rdms_searchlight.m in the rdms/ directory
%

clear all;
%close all;

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
%dirname = 'rdms/control';
%dirname = 'rdms/M1M2M1_4mm';
%dirname = 'rdms/M1M2M1_4mm_nosmooth';
dirname = 'rdms/M1M2M1_4mm_ordered';

EXPT = context_expt();
rsa_idx = 1;
model_idx = 1;



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
tmap2 = tmap;
pmap = tmap;
V.fname = fullfile(dirname, ['searchlight_tmap_', model, '_', event, '.nii']); % change immediately!

maxt = 0;



rsadir = fullfile(EXPT.rsadir,['rsa',num2str(rsa_idx)]);

% create sample rsa
rsa = EXPT.create_rsa(rsa_idx, 1);

files_ccnl = {};
% compute t-map, if it hasn't been computed already
%if ~exist(V.fname, 'file')
    df = NaN;
    files = dir(rsadir);
    for i = 1:length(files)
        file = files(i).name;
        if startsWith(file, 'searchlight_') && endsWith(file, '.mat')
            files_ccnl = [files_ccnl {file}];
            %df = size(all_subject_rhos, 3) - 1;
        end
    end

    %V.descrip = sprintf('SPM{T_[%d.0]}', df); % hack to write the degrees of freedom, to allow thresholding in bspmview

    % save tmap
    %V.fname
    %spm_write_vol(V, tmap);
%end



%% load all the searchlight results
%
files = dir(dirname);
ccnl_idx = 0;
for i = 1:length(files)
    file = files(i).name;
    if startsWith(file, 'searchlight_') && endsWith(file, '.mat')
        disp(['Loading ', file]);
        table_T = [];
        load(fullfile(dirname, file));

        if strcmp(file, 'searchlight_weights_217001-217783.mat')
            continue % hack; this one is absent from new rsa TODO fix 
        end

        ccnl_idx = ccnl_idx + 1;
        disp(['Loading also ', files_ccnl{ccnl_idx}]);
        load(fullfile(rsadir, files_ccnl{ccnl_idx}));


        % when using _nosmooth, sometimes not all searchlights are good
        % => have to keep 
        %
        for j = 1:length(x)
            if ~isequal(event, events{j}), continue; end
            row_idx = j;
            col_idx = find(ismember(model_names, model));
            if col_idx <= size(table_T, 2) % in case we're using a model that hasn't been computed in all batches, e.g. weightsPosterior
                tmap(x(j), y(j), z(j)) = table_T(row_idx, col_idx);

                tmap2(cor(j,1), cor(j,2), cor(j,3)) = T(row_idx, 1);

                assert(col_idx == 1);
                assert(immse(T(row_idx, 1), table_T(row_idx, col_idx)) < 1e-10);
            end
        end
        clear events; % for next iteration
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

%spm_write_vol(V, tmap);
spm_write_vol(V, tmap2);

%% visualize
%
EXPT = context_expt();
struc = fullfile('masks','mean.nii');

disp(V.fname)

bspmview(V.fname, struc);
