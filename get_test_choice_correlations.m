

% correlate model choices with subject choices on the test trials
%
function [r, p] = get_test_choice_correlations(params_file, params_idx, which_structures)
    utils; % include some nifty lambdas

    [data, metadata, simulated] = simulate_subjects_helper(true, params_file, params_idx, which_structures);

    %
    % Choice probabilities in test phase for SUBJECTS
    %

    Ms = [];
    SEMs = [];
    for context = metadata.contextRoles
        which = data.which_rows & data.isTrain == 0 & strcmp(data.contextRole, context);

        x1c1 = strcmp(data.response.keys(which & data.cueId == 0 & data.contextId == 0), 'left');
        x1c2 = strcmp(data.response.keys(which & data.cueId == 0 & data.contextId == 2), 'left');
        x2c1 = strcmp(data.response.keys(which & data.cueId == 2 & data.contextId == 0), 'left');
        x2c2 = strcmp(data.response.keys(which & data.cueId == 2 & data.contextId == 2), 'left');

    %    M = mean([x1c1 x1c2 x2c1 x2c2]);
    %    SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
        M = get_means(x1c1, x1c2, x2c1, x2c2);
        SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
        Ms = [Ms; M];
        SEMs = [SEMs; SEM];
    end

    subject_Ms = Ms; % for stats

    %
    % TRUE Choice probabilities in test phase for MODEL
    %

    Ms = [];
    SEMs = [];
    for context = metadata.contextRoles
        which = data.which_rows & data.isTrain == 0 & strcmp(data.contextRole, context);

        x1c1 = simulated.pred(which & data.cueId == 0 & data.contextId == 0);
        x1c2 = simulated.pred(which & data.cueId == 0 & data.contextId == 2);
        x2c1 = simulated.pred(which & data.cueId == 2 & data.contextId == 0);
        x2c2 = simulated.pred(which & data.cueId == 2 & data.contextId == 2);

        %M = mean([x1c1 x1c2 x2c1 x2c2]);
        %SEM = std([x1c1 x1c2 x2c1 x2c2]) / sqrt(length(x1c1));
        M = get_means(x1c1, x1c2, x2c1, x2c2);
        SEM = get_sems(x1c1, x1c2, x2c1, x2c2);
        Ms = [Ms; M];
        SEMs = [SEMs; SEM];
    end

    model_Ms = Ms; % for stats

    % correlate average subject choices with model choices 
    %
    [r, p] = corrcoef(subject_Ms(:), model_Ms(:));
    r = r(1,2);
    p = p(1,2);
    
    figure;
    scatter(subject_Ms(:)+rand(size(subject_Ms(:)))*0.01, model_Ms(:)+rand(size(subject_Ms(:)))*0.01);
    xlabel('human');
    ylabel('model');
    lsline;
    save('wtf.mat');
end
