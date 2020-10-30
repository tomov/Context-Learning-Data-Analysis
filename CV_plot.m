% plot results of CV

glmodels = which_structuress;
glmodels{1} = 'ideal';

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

m = 3;
for g = 1:3
    [h,p,ci,stat] = ttest(atanh(r_avg{m}(:,4)), atanh(r_avg{m}(:,g)))
end
