% plot results of CV

glmodels = which_structuress;
glmodels{1} = 'ideal';

figure;

for m = 1:length(masks)

    x = [];
    y = [];
    for g = 1:length(glmodels)
        x = [x g*ones(1, length(goodSubjs))];
        %y = [y mse{m}(:,g)'];
        %y = [y atanh(r_avg{m}(:,g))'];
        y = [y avg_r{m}(:,g)'];
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
    %[h,p,ci,stat] = ttest(atanh(r_avg{m}(:,4)), atanh(r_avg{m}(:,g)));
    [h,p,ci,stat] = ttest(atanh(avg_r{m}(:,4)), atanh(avg_r{m}(:,g)));
    %[h,p,ci,stat] = ttest(mse{m}(:,4), mse{m}(:,g));
    fprintf('%s vs. %s: t(%d) = %.4f, p = %.4f\n', glmodels{4}, glmodels{g}, stat.df, stat.tstat, p);
    %[p,h,stat] = signrank(mse{m}(:,4), mse{m}(:,g));
    %fprintf('%s vs. %s: p = %.4f\n', glmodels{4}, glmodels{g}, p);
end

for g = 1:4
    %[h,p,ci,stat] = ttest(atanh(r_avg{m}(:,g)));
    [h,p,ci,stat] = ttest(atanh(avg_r{m}(:,g)));
    %[h,p,ci,stat] = ttest(mse{m}(:,g));
    fprintf('%s vs. 0: t(%d) = %.4f, p = %.4f\n', glmodels{g}, stat.df, stat.tstat, p);
    %[p,h,stat] = signrank(avg_r{m}(:,g));
    %[p,h,stat] = signrank(mse{m}(:,g));
    %fprintf('%s vs. 0: p = %.4f\n', glmodels{g}, p);
end
