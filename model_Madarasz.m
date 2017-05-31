% Script that simulates Madarasz et al. (2016) Nature Neuroscience
%
close all;
clear all;

params = [1 1];
which_structures = [1 1 1 0];
struct_names = {'M1 (irr)', 'M2 (mod)', 'M3 (add)'};

conditions = {'control-1', 'control-2', 'intermixed', 'pairings-first'};

x{1} = [0; 0; 1; 0; 0; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0]; % CS
r{1} = [0; 0; 1; 0; 0; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0]; % US
%k{1} = ones(size(x{1})); % context A
k{1} = [1; 1; 1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2; 2; 2];

x{2} = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % CS
r{2} = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % US
%k{2} = ones(size(x{1})); % context A
k{2} = [1; 1; 1; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2];

x{3} = [0; 0; 1; 0; 0; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0]; % CS
r{3} = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]; % US
%k{3} = ones(size(x{1})); % context A
k{3} = [1; 1; 1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2; 2; 2];

x{4} = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % CS
r{4} = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1]; % US
%k{4} = ones(size(x{1})); % context A
k{4} = [1; 1; 1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2; 2; 2];

test_x = [1; 0];
test_k = [3; 1];

P_ns = {};
ww_ns = {};
Ps = {};
wws = {};

Ms = nan(numel(test_k), numel(conditions));
for cond_id = 1:numel(conditions)
    condition = conditions{cond_id};

    [choices, P_n, ww_n, P, ww] = model_train_Madarasz(x{cond_id}, k{cond_id}, r{cond_id}, params, which_structures, false);
    [test_choices] = model_test_Madarasz(test_x, test_k, P_n, ww_n, params);
    
    Ms(1, cond_id) = test_choices(1);
    Ms(2, cond_id) = test_choices(2);
    
    P_ns{cond_id} = P_n;
    ww_ns{cond_id} = ww_n;
    Ps{cond_id} = P;
    wws{cond_id} = ww;
end

figure;
barweb(Ms, zeros(size(Ms)), 1, {'Tone memory', 'Context memory'}, 'CR (freezing)');
legend(conditions);

figure;
P_ns = cell2mat(P_ns');
P_ns = P_ns(:, logical(which_structures));
barweb(P_ns, zeros(size(P_ns)), 1, conditions, 'Posterior over structures');
legend(struct_names);