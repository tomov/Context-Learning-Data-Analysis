% Helper f'n that predicts the value on the current trial
%
function [V, vals] = model_value_helper(x, w, stimulus, context, P)
    for i = 1:numel(P)
        if i == 2 
            % for M2, use weights for current context
            %
            vals(i) = x{i}' * w{i}(:,context);
        elseif i == 5
            % for M5 = M2', use weights for current cue
            %
            vals(i) = x{i}' * w{i}(:,find(stimulus));
        else
            vals(i) = x{i}' * w{i};
        end
    end
    V = sum(vals .* P);
end
