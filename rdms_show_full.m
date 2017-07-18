function rdms_show_full(rows, cols, table_Rho, table_P, all_vs_all)

% Visualize full second-order RDM
% Show full second-order RSA matrix with correlation coefficients + another one with the p-values
%
% INPUT:
% rows, cols = struct arrays of RDMs as input to rdms_second_order()
% table_Rho, table_P = correlation coefficients and p-values as output by rdms_second_order()
% all_vs_all = whether we're comparing all RDMs vs. themselves, or Neural vs. Model RDMs
%


tabs = {table_Rho, table_P};
titles = {'Representational similarity match (Spearman rank correlation)', 'P-value (one-sample t-test)'};

for i = 1:numel(tabs)
    figure;
    if i == 2
        t = tabs{i};
        if all_vs_all
            imagesc(log10(t), [-16 0.1]);
        else
            imagesc(log10(t));
        end    
        c = colorbar;
        yt = get(c, 'Ytick');
        set(c, 'YTickLabel', arrayfun(@(x) ['10\^', num2str(x)], yt, 'UniformOutput', false));
    else
        if all_vs_all
            imagesc(tabs{i}, [0 0.4]);
        else
            imagesc(tabs{i});
        end
        %cols=colorScale([0 0.5 1; 0.5 0.5 0.5; 1 0 0],256);
        %colormap(cols); colorbar;
        colorbar;
    end

    xticklabels({cols.name});
    set(gca, 'xtick', 1:numel({cols.name}));
    xtickangle(60);
    xlabel('Neural model');

    yticklabels({rows.name});
    set(gca, 'ytick', 1:numel({rows.name}));
    ylabel('ROI_{event}, t = trial onset, f = feedback onset');
    
    title(titles{i});
end

