% Play around with linear mixed effects (LME) modeling
%

%% Set up some data for two "subjects"
%
x1 = 1:10;
x2 = 0.5:10.5;
y = [];
y = [y; glmval([6.8 1.5]', x1', 'identity')]; % subj 1
y = [y; glmval([-0.5 2.9]', x2', 'identity')]; % subj 2
x = [x1 x2]';
y = y + rand(size(y)) * 3; % some random noise
id = [ones(numel(x1), 1); 2 * ones(numel(x2), 1)];

tbl = array2table([y x id], 'VariableNames', {'Y', 'X', 'subj'});
disp(tbl);



%% Run LME but with fixed effects only
%
lme = fitlme(tbl, 'Y ~ X');

[beta, betanames, stats] = fixedEffects(lme);
res = residuals(lme);

% Plot fit + residuals
%
scatter(x, y);
hold on;
z = glmval(beta, x, 'identity');
plot(x, z);
for i = 1:numel(x)
    plot([x(i) x(i)], [z(i) z(i) + res(i)], 'Color', 'green');
end

legend({'data', 'fixed effects', 'residuals'});
hold off;



%% Now run LME with random effects
%
lme = fitlme(tbl, 'Y ~ X + (X|subj)');

[beta, betanames, stats] = fixedEffects(lme);
res = residuals(lme);

[rbeta, rbetanames, rstats] = randomEffects(lme);

% Plot fit + residuals
%
scatter(x, y);
hold on;
z = glmval(beta, x, 'identity');
z1 = glmval(beta + rbeta(1:2), x1, 'identity');
z2 = glmval(beta + rbeta(3:4), x2, 'identity');

% plot fixed effect
plot(x, z);

% plot random effect
plot(x1, z1);
plot(x2, z2);

% plot residuals
z = [z1; z2];
for i = 1:numel(x)
    plot([x(i) x(i)], [z(i) z(i) + res(i)], 'Color', 'green');
end

legend({'data', 'fixed effects', 'random effects, subj 1', 'random effects, subj 2', 'residuals'});
hold off;
