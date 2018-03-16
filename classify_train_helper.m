function [classifier, outputs, accuracy] = classify_train_helper(method, inputs, targets, runs, trials, subjs, outFilename)

% helper for classify_train.m
%

accuracy = NaN; % TODO make all of them output it


tic
disp('training classifier...');

%
% Fit them
%

switch method
    
    case 'patternnet' % neural network classifier
        
        % patternnet wants column feature vectors. I.e. each data point is a column
        % so we have to rotate it ...
        %
        inputs = inputs'; % ugh MATLAB
        targets = targets';

        % from https://github.com/tomov/food-recognition/blob/master/neural_train.m

        % Create a Pattern Recognition Network
        hiddenLayerSize = 8; % TODO param
        net = patternnet(hiddenLayerSize);

        % Set up Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;
        net.trainParam.showWindow = false; % don't show GUI on NCF

        % Train the Network
        [net,tr] = train(net,inputs,targets);

        % Test the Network
        outputs = net(inputs);
        errors = gsubtract(targets,outputs);
        performance = perform(net,targets,outputs);

        % View the Network
        %view(net)

        % View confusion matrix
        [c,cm,ind,per] = confusion(targets,outputs);

        % patternnet wants column feature vectors, so we rotated those
        % but the rest of the code expects them to be rotated
        % so we undo that...
        % 
        inputs = inputs';
        targets = targets';
        outputs = outputs';

        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate = %.2f%%\n', accuracy);    

        classifier = net;

        % TODO these break for some reason...
        %
        %{
        A = cm';
        B = bsxfun(@rdivide,A,sum(A));
        C = imresize(1 - B, 15, 'method', 'box');
        imshow(C);
        xlabel('Targets');
        ylabel('Outputs');


        [b, ib] = sort(diag(B));
        sprintf('%.2f %% bottom 10', mean(b(1:10) * 100))
        sprintf('%.2f %% top 10', mean(b(41:50) * 100))
        sprintf('%.2f %% correct', (1 - c) * 100)
        %}
    
        
    case 'glmnet' % multinomial GLM classifier
    
        opts.alpha = 1; % 0 = ridge penalty; 1 = lasso penalty (force betas to zero); default = 1
        opts.mtype = 'ungrouped'; % 'grouped' = for multinomial, all betas are in our out together; default = 'ungrouped'
        opts.nlambda = 1000; % # lambdas to use; default = 100
        opts.lambda_min = 0.00000001; % as a fraction of lambda_max which is derived from the data; default = 0.0001
        options = glmnetSet(opts);

        % x = inputs
        % y = targets
        %
        fitObj = glmnet(inputs, targets, 'multinomial', options);
        glmnetPrint(fitObj);

        % no need to do that manually here -- cvglmnet does it
        % automatically
        %
        %[mses, msesems] = glmnetKFoldCrossValidation(inputs, targets, fitObj, 'multinomial', 'response', 4);
        %[~, lambda_idx] = min(mses); % pick lambda with smallest MSE
        %lambda = fitObj.lambda(lambda_idx);

        outputss = glmnetPredict(fitObj, inputs, fitObj.lambda, 'response');

        accuracies = nan(size(outputss, 3), 1);
        for l = 1:size(outputss, 3) % for all lambdas
            outputs = outputss(:,:,l);
            accuracy = classify_get_accuracy(outputs, targets);
            accuracies(l) = accuracy;
            fprintf('Success rate for %d (lambda = %.4f) = %.2f%%\n', l, fitObj.lambda(l), accuracy);
        end

        classifier = fitObj;
   

        
    case 'cvglmnet' % cross-validated multinomial GLM classifier
    
        opts.alpha = 1; % 0 = ridge penalty; 1 = lasso penalty (force betas to zero); default = 1
        opts.mtype = 'ungrouped'; % 'grouped' = for multinomial, all betas are in our out together; default = 'ungrouped'
        opts.nlambda = 1000; % # lambdas to use; default = 100
        opts.lambda_min = 0.00000001; % as a fraction of lambda_max which is derived from the data; default = 0.0001
        options = glmnetSet(opts);

        nfolds = 10;
        c = cvpartition_runs(nfolds, subjs, runs, trials, inputs);
        foldid = c.Impl.indices;

        % each run is a separate fold
        %
        %{
        if isempty(foldid)
            [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
            foldid = data.runId(which_rows);
            foldid = arrayfun(@(x) find(x == runs), foldid);
        end
        %}
        % each (subject, run) is a separate fold 
        %
        %{
        [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
        cur_participant = '';
        cur_run = NaN;
        cur_fold = 0;
        foldid = [];
        for i = find(which_rows)'
            if ~isequal(data.participant{i}, cur_participant) || ~isequal(data.runId(i), cur_run)
                cur_fold = cur_fold + 1;
            end
            foldid = [foldid; cur_fold];
            cur_participant = data.participant{i};
            cur_run = data.runId(i);
        end
        %}
        assert(length(foldid) == size(inputs, 1));
        disp('folds:');
        disp(foldid);

        
        % x = inputs
        % y = targets
        %
        parallel = false;
        keep = true;
        CVfit = cvglmnet(inputs, targets, 'multinomial', options, 'deviance', [], foldid, parallel, keep);
        disp(CVfit);

        outputs = cvglmnetPredict(CVfit, inputs, CVfit.lambda_1se, 'response');

        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate (lambda = %.4f) is %.2f%%\n', CVfit.lambda_1se, accuracy);

        classifier = CVfit;


    case 'fitcecoc' % multinomial SVM classifier, cross-validated

        kfolds = 10;
        c = cvpartition_runs(kfolds, subjs, runs, trials, inputs);

        [~, labels] = max(targets, [], 2);
        Mdl = fitcecoc(inputs, labels, 'CVPartition', c, 'FitPosterior', 1, 'Verbose', 2);

        [out_labels, negloss, cost, outputs] = kfoldPredict(Mdl);

        %accuracy = mean(outputs == labels);  lame -- don't use this
        %accuracy = mean(sum(outputs .* targets, 2));
        accuracy = classify_get_accuracy(outputs, targets);
        fprintf('Success rate is %.2f%%\n', accuracy);

        classifier = Mdl;


    case 'mnrfit'
        % 
        %
        
    otherwise
        assert(false, 'should be one of the above');
end

toc

% Save everything except for inputs (too big)
%
fprintf('SAVING to %s\n', outFilename);
save(outFilename,'-regexp','^(?!(inputs)$).');

end



% ensure folds contain whole runs, i.e. no run is split across training & test
% this is to ensure cross-run predictions
% b/c of temporal autocorrelation within runs, it is very easy to classify conditions within runs
%
function c = cvpartition_runs(k, subjs, runs, trials, inputs)
    c_runs = cvpartition(numel(subjs) * numel(runs), 'Kfold', k);
    assert(size(inputs, 1) == numel(subjs) * numel(runs) * numel(trials));
    c = cvpartition(size(inputs, 1), 'Kfold', k);
    i = repmat(c_runs.Impl.indices, 1, numel(trials));
    i = i'; 
    i = i(:);
    c.Impl.indices = i;
    c.Impl.TestSize = accumarray(i, 1)';
    c.Impl.TrainSize = size(i, 1) - c.Impl.TestSize;
end
