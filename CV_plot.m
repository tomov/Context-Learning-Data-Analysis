% plot results of CV

figure;

for m = 1:length(masks)

    x = [];
    y = [];
    for g = 1:length(glmodels)
        x = [x g*ones(1, length(goodSubjs))];
        y = [y mse{m}(:,g)'];
        %y = [y atanh(r_avg{m}(:,g))'];
    end 

    subplot(length(masks), 1, m);
    swarmchart(x, y);
    title(masks{m}, 'interpreter', 'none');
    ylabel('MSE');
    xticks(1:length(glmodels));
    xticklabels(glmodels);
    xlabel('GLM');
end 
