function extract_results_table(varargin)
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
% V = SPM volume of the t-map, with the filename changed so we don't
%     overwrite it accidentally
% Y = the actual t-map
% C = volume with cluster size for each voxel
% CI = volume with cluster index for each voxel <-- that's the name of the
%      game; 
% region = labels for the peak voxels
% extent = size of the cluster
% stat = t-statistic for the peak voxels
% mni = MNI coordinates of peak voxels
% cor = coordinates of peak voxel in native space (can be used as indices
%       in the C and CI volumes)
% results_table = what Save Results Table in bspmview would spit out 
% 


atlas_name = varargin{1};
method = varargin{2};
assert(ismember(atlas_name, {'AnatomyToolbox', 'AAL2', 'HarvardOxford-maxprob-thr0', 'Talairach', 'Brodmann'}));
assert(ismember(method, {'peak', 'vote', 'all'}));

[V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = extract_clusters(varargin{3:end});

atlas_dirpath = 'atlases';


[atlaslabels, atlas] = bspm_setatlas(spmT, atlas_dirpath, atlas_name);
atlas = reshape(atlas, size(Y)); % voxel --> num
map = containers.Map(atlaslabels.id, atlaslabels.label); % num --> label

new_region = {};

for i = 1:size(region, 1) 
    % get peak voxel
    x = cor(i,1);
    y = cor(i,2);
    z = cor(i,3);
    assert(immse(stat(i), Y(x,y,z)) < 1e-6);

    % get cluser mask and voxels
    clust_idx = CI(x,y,z);
    mask = CI == clust_idx;
    voxels = find(mask);
    assert(numel(voxels) == extent(i)); % <-- doesn't work with +/- ; do one at a time

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
                new_region{i} = strjoin(rois, ' \\\\\n & ');
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
    
    fprintf('%s & %s & %s & %d & %.3f & %d %d %d \\\\\n', sign, new_region{i}, '??', extent(i), stat(i), mni(i,1), mni(i,2), mni(i,3));
    if isequal(method, 'all')
        fprintf('\\hline\n');
    end

    
end

save('shit.mat');

%new_region'
