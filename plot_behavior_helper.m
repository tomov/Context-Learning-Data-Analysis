% helper f'n to plot behavioral results
%
function plot_behavior_helper(data, metadata, simulated)

    % Choice probabilities for model
    %
    model_means = [];
    model_sems = [];
    for condition = metadata.contextRoles
        which_rows = data.which_rows & ~data.isTrain & strcmp(data.contextRole, condition);
        
        x1c1 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 0);
        x1c3 = simulated.pred(which_rows & data.cueId == 0 & data.contextId == 2);
        x3c1 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 0);
        x3c3 = simulated.pred(which_rows & data.cueId == 2 & data.contextId == 2);

        model_means = [model_means; mean(x1c1) mean(x1c3) mean(x3c1) mean(x3c3)];
        model_sems = [model_sems; std(x1c1)/sqrt(numel(x1c1)) std(x1c3)/sqrt(numel(x1c3)) std(x3c1)/sqrt(numel(x3c1)) std(x3c3)/sqrt(numel(x3c3))];
    end

    model_means

    % Plot model choices probabilities
    %
    h = bar(model_means);

    hold on;

    if nargin == 2 & ~isempty(model_sems)
        xs = sort([h(1).XData + h(1).XOffset, ...
                   h(2).XData + h(2).XOffset, ...
                   h(3).XData + h(3).XOffset, ...
                   h(4).XData + h(4).XOffset]);
        model_means = model_means'; model_means = model_means(:);
        model_sems = model_sems'; model_sems = model_sems(:);        
        errorbar(xs, model_means, model_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0], 'LineWidth', 1, 'Color', [0 0 0], 'AlignVertexCenters', 'off');
    end

    % plot both pilot and fmri subjects' choices
    %
    group_x_offs = [-0.02, 0.02];
    group_color = {[0.5 0.5 0.5], [0 0 0]};
    for group = 1:2

        if group == 1
            % pilot subjects
            [data, metadata] = load_data(fullfile('data', 'pilot.csv'), false);
            disp('\n\n\n --------------------------------------- PILOT -----------------\n\n');
        else
            % fmri subjects
            [data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
            disp('\n\n\n --------------------------------------- fMRI -----------------\n\n');
        end
        
        % Choice probabilities for human subjects
        %
        human_choices = [];        
        for condition = metadata.contextRoles
            which_rows = data.which_rows & ~data.isTrain & strcmp(data.contextRole, condition);

            x1c1 = []; x1c3 = []; x3c1 = []; x3c3 = [];
            for who = metadata.subjects
                x1c1 = [x1c1, data.chose_sick(which_rows & data.cueId == 0 & data.contextId == 0 & strcmp(data.participant, who))];
                x1c3 = [x1c3, data.chose_sick(which_rows & data.cueId == 0 & data.contextId == 2 & strcmp(data.participant, who))];
                x3c1 = [x3c1, data.chose_sick(which_rows & data.cueId == 2 & data.contextId == 0 & strcmp(data.participant, who))];
                x3c3 = [x3c3, data.chose_sick(which_rows & data.cueId == 2 & data.contextId == 2 & strcmp(data.participant, who))];
            end
            assert(isequal(size(x1c1), [metadata.runsPerContext metadata.N]));
            assert(isequal(size(x1c3), [metadata.runsPerContext metadata.N]));
            assert(isequal(size(x3c1), [metadata.runsPerContext metadata.N]));
            assert(isequal(size(x3c3), [metadata.runsPerContext metadata.N]));

            human_choices = [human_choices; mean(x1c1); mean(x1c3); mean(x3c1); mean(x3c3)];
        end

        [human_sems, human_means] = wse(human_choices');


        % Get x-coordinates of bar centers and "flatten" the human choice
        % data so we can overlay it nicely
        %
        xs = sort([h(1).XData + h(1).XOffset, ...
                   h(2).XData + h(2).XOffset, ...
                   h(3).XData + h(3).XOffset, ...
                   h(4).XData + h(4).XOffset]);
        human_means = human_means'; human_means = human_means(:);
        human_sems = human_sems'; human_sems = human_sems(:);        

        % Plot human choices
        %
        errorbar(xs + group_x_offs(group), human_means, human_sems, '.', 'MarkerSize', 1, 'MarkerFaceColor', group_color{group}, 'LineWidth', 1, 'Color', group_color{group}, 'AlignVertexCenters', 'off');
        % the markers in the errorbar plot are misaligned from the error
        % bars -- this hack adjusts them
        xs_adjust = xs + 0.003;
        %{
        xs_adjust(1) = xs_adjust(1) + 0.0035;
        xs_adjust(2) = xs_adjust(2) + 0.003;
        xs_adjust(3) = xs_adjust(3) + 0.003;
        xs_adjust(4) = xs_adjust(4) + 0.0025;
        xs_adjust(5) = xs_adjust(5) + 0.004;
        xs_adjust(6) = xs_adjust(6) + 0.004;
        xs_adjust(7) = xs_adjust(7) + 0.0025;
        xs_adjust(8) = xs_adjust(8) + 0.0025;
        xs_adjust(9) = xs_adjust(9) + 0.003;
        xs_adjust(10) = xs_adjust(10) + 0.004;
        xs_adjust(11) = xs_adjust(11) + 0.003;
        xs_adjust(12) = xs_adjust(12) + 0.0025;
        %}
        plot(xs_adjust + group_x_offs(group), human_means, 'o', 'MarkerSize', 5, 'MarkerFaceColor', group_color{group}, 'Color', group_color{group});
    end
    
    
    
    hold off;
    
    %xticklabels({'Irrelevant training', 'Modulatory training', 'Additive training'});
    xticklabels({'Irrelevant context', 'Modulatory context', 'Irrelevant cue'});
    ylabel('Choice probability');
    legend({'x_1c_1', 'x_1c_3', 'x_3c_1', 'x_3c_3'}, 'Position', [0.07 -0.095 1 1]);
    ylim([0 1.1]);
    set(gca,'fontsize',13);
    
    %print(gcf, 'untitled.pdf', '-dpdf', '-bestfit');
    
    text(0.1, 1.25, 'B', 'FontSize', 20, 'FontWeight', 'bold')
end

