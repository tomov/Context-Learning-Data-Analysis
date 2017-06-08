% RDMs
%

% Load the subject behavioral data
%
[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());

rois = {'hippocampus', 'ofc', 'vmpfc', 'striatum', 'pallidum', 'bg', ...
        'v1', 'm1', 's1', 'fusiform', 'angular', 'mid_front', 'dl_sup_front'};

masks = cellfun(@(x) fullfile('masks', [x, '.nii']), rois, 'UniformOutput', false);
distance_measures = {'correlation', 'cosine', 'euclidean'};
regressor_prefixes = {'trial_onset', 'feedback_onset'};

for mask = masks
    mask = mask{1}
    
    for distance_measure = distance_measures
        distance_measure = distance_measure{1}
        
        for regressor_prefix = regressor_prefixes
            regressor_prefix = regressor_prefix{1}

            % Load the RDMs for that mask
            %
            [subjectRDMs, avgSubjectRDM, concatSubjectRDMs] = get_neural_rdms(mask, distance_measure, regressor_prefix, data, metadata);

            % Display the RDMs
            %
            %showRDMs(concatSubjectRDMs, 2);
            %showRDMs(avgSubjectRDM, 1);
        end
    end
end

