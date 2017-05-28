% A bunch of convenient lambdas
%

sem = @(x) std(x) / sqrt(length(x));

softmax = @(V_n, inv_temp) 1 ./ (1 + exp(-2 * inv_temp * V_n + inv_temp));

% b/c sometimes they're vectors of size 1 == scalars, so can't do mean([a b c d e]) 
%
get_means = @(x1c1, x1c2, x2c1, x2c2) [mean(x1c1) mean(x1c2) mean(x2c1) mean(x2c2)];
get_sems = @(x1c1, x1c2, x2c1, x2c2) [std(x1c1) / sqrt(length(x1c1)) std(x1c2) / sqrt(length(x1c2)) std(x2c1) / sqrt(length(x2c1)) std(x2c2) / sqrt(length(x2c2))];

