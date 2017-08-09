function table = extract_results_table(varargin)
% Wrapper around extract_clusters that also uses different atlases
%
% Given a contrast, extract all the activation clusters from the t-map after cluster FWE
% correction.
% Then it uses different atlases to figure out the names of the different activation clusters.
% As a bonus, prints out the table in LaTeX format.
%
% INPUT:
% atlas_name = which atlas to use, e.g. 'AAL2' or 'AnatomyToolbox'
% method = 'peak', 'vote', or 'all'
% the rest of the parameters are the same as extract_clusters()
%
% OUTPUT:
% table = a nice table with all the stuffs
% 
% EXAMPLE:
% table = extract_results_table('AnatomyToolbox', 'peak', context_expt(), 154, 'KL_structures', 0.001, '+', 0.05);
% table = extract_results_table('Brodmann', 'vote', context_expt(), 154, 'KL_structures', 0.001, '+', 0.05);

atlas_dirpath = 'atlases';

atlas_name = varargin{1};
method = varargin{2};
assert(ismember(atlas_name, {'AnatomyToolbox', 'AAL2', 'HarvardOxford-maxprob-thr0', 'Talairach', 'Brodmann'}));
assert(ismember(method, {'peak', 'vote', 'all'}));

% extract the clusters
%
[V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = extract_clusters(varargin{3:end});

% load atlas
%
[atlaslabels, atlas] = bspm_setatlas(spmT, atlas_dirpath, atlas_name);
atlas = reshape(atlas, size(Y)); % voxel --> num
map = containers.Map(atlaslabels.id, atlaslabels.label); % num --> label

% load BA atlas separately
%
[~, BAatlas] = bspm_setatlas(spmT, atlas_dirpath, 'Brodmann');
BAatlas = reshape(BAatlas, size(Y)); % voxel --> BA

% label the clusters
%
new_region = {};

for i = 1:size(region, 1) 
    % get peak voxel
    x = cor(i,1);
    y = cor(i,2);
    z = cor(i,3);
    assert(immse(stat(i), Y(x,y,z)) < 1e-6);
    
    % set brodmann area
    if BAatlas(x, y, z)
        BA = num2str(BAatlas(x, y, z));
    else
        BA = '';
    end

    % get cluser mask and voxels
    clust_idx = CI(x,y,z);
    mask = CI == clust_idx;
    voxels = find(mask);
    if ~isequal(method, 'peak')
        assert(numel(voxels) == extent(i)); % <-- doesn't work with +/- ; do one at a time
    end

    % see which region each voxel "votes" for
    %
    [x, y, z] = ind2sub(size(mask), voxels);
    votes = zeros(1, max(atlas(:)));
    sign_votes = zeros(1, max(atlas(:))); % whether it's L or R on average
    for j = 1:numel(x)
        if atlas(x(j),y(j),z(j)) > 0
            idx = atlas(x(j),y(j),z(j));
            votes(idx) = votes(idx) + 1;
            
            m = cor2mni([x(j) y(j) z(j)], V.mat);
            if m(1) < 0
                sign_votes(idx) = sign_votes(idx) - 1; % L voxel
            else
                sign_votes(idx) = sign_votes(idx) + 1; % R voxel
            end
        end
    end

    % get the cluster label(s)
    %
    switch method
        case 'peak'
            % just use the peak voxel
            %
            if isKey(map, atlas(cor(i,1),cor(i,2),cor(i,3)))
                label = map(atlas(cor(i,1),cor(i,2),cor(i,3)));
                new_region{i} = atlas_label_to_roi_name(atlas_name, label, mni(i,:)); % convert to nice anatomical name
            else
                new_region{i} = '';
            end
            
        case 'vote'
            % get the most popular anatomical region in the cluster
            % each voxel gets one vote
            %
            if sum(votes) == 0
                new_region{i} = ''; % no region :(
            else
                [~, idx] = max(votes);
                label = map(idx);
                new_region{i} = atlas_label_to_roi_name(atlas_name, label, [sign_votes(idx) 0 0]); % convert to nice anatomical name
            end
            
        case 'all'
            % List the regions of all voxels
            %
            if sum(votes) == 0
                new_region{i} = ''; % no region :(
            else
                [~, idxs] = find(votes);
                [~, sorted_idxs] = sort(votes(idxs), 'descend');

                rois = {};
                for j = 1:numel(sorted_idxs)
                    idx = idxs(sorted_idxs(j));
                    assert(votes(idx) > 0);

                    label = map(idx);
                    roi = atlas_label_to_roi_name(atlas_name, label, [sign_votes(idx) 0 0]);
                    roi = sprintf('%s (%.2f\\%%)', roi, 100 * votes(idx) / sum(votes));
                    rois = [rois; {roi}];
                end
                new_region{i} = rois;
            end
            
        otherwise
            assert(false, 'method should be one of these');
    end
    
    % print for LaTeX
    %
    sign = '';
    if i == 1 || (stat(i) < 0) ~= (stat(i - 1) < 0) % sign changed
        if stat(i) < 0
            sign = 'Negative';
        else
            sign = 'Positive';
        end
    end
    
    region_latex = new_region{i};
    if iscell(new_region{i}) % multiple regions in this cluster
        region_latex = strjoin(new_region{i}, ' \\\\\n & ');
    else
        region_latex = new_region{i};
    end
    fprintf('%s & %s & %s & %d & %.3f & %d %d %d \\\\\n', sign, region_latex, BA, extent(i), stat(i), mni(i,1), mni(i,2), mni(i,3));
    if isequal(method, 'all')
        fprintf('\\hline\n');
    end
    
    % output table
    %
    table{i, 1} = sign;
    table{i, 2} = new_region{i};
    table{i, 3} = BA;
    table{i, 4} = extent(i);
    table{i, 5} = stat(i);
    table{i, 6} = mni(i,1);
    table{i, 7} = mni(i,2);
    table{i, 8} = mni(i,3);

    
end

save('shit.mat');

%new_region'
