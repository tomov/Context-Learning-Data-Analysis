function [V, Y, C, CI, region, extent, stat, mni, cor, results_table] = extract_results_table(varargin)
% Wrapper around bspm_extract_clusters
% Given a contrast, extract all the activation clusters from the t-map after cluster FWE
% correction. Uses the same logic as bspmview,
% and as a sanity check prints the results table -- should be the same as
% the one from bspmview.
% Then it uses different atlases to figure out the names of the different activation clusters.
% As a bonus, prints out the table in LaTeX format.
%
% INPUT:
% same as extract_clusters()
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

method = 'peak';

[V, Y, C, CI, region, extent, stat, mni, cor, results_table, spmT] = extract_clusters(varargin{:});

atlas_dirpath = '/Users/momchil/Dropbox/Research/libs/bspmview/supportfiles';

atlas_name = 'HarvardOxford-maxprob-thr0';
atlas_name = 'AAL2';
atlas_name = 'AnatomyToolbox';

%atlas_names = {'AnatomyToolbox', 'AAL2', 'HarvardOxford-maxprob-thr0'}

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

    % get cluser voxels
    clust_idx = CI(x,y,z);
    mask = CI == clust_idx;

    voxels = find(mask);
    assert(numel(voxels) == extent(i)); % <-- doesn't work with +/- ; do one at a time

    % get peak voxel atlas label 
    if isKey(map, atlas(x,y,z))
        new_region{i} = map(atlas(x, y, z));
    else
        new_region{i} = '';
    end
    
    % Convert labels to nice anatomical regions
    %
    switch atlas_name
        case 'AAL2'
            new_region{i} = aal2_label_to_roi_name(new_region{i}, mni(i,:));
            
        case 'AnatomyToolbox'
            hemi = '';
            if startsWith(new_region{i}, 'L ')
                new_region{i} = new_region{i}(3:end);
                hemi = 'L';
            elseif startsWith(new_region{i}, 'R ')
                new_region{i} = new_region{i}(3:end);
                hemi = 'R';
            end
            
            space = find(new_region{i} == ' ' | new_region{i} == '-');
            if ~isempty(space)
                new_region{i} = [new_region{i}(1:space), lower(new_region{i}(space+1:end))];
            end
            
            if ~isempty(hemi)
                new_region{i} = [new_region{i}, ' (', hemi, ')'];
            end
            
        otherwise
            assert(false, 'Should be one of the above');
    end

    sign = '';
    if i == 1 || (stat(i) < 0) ~= (stat(i - 1) < 0) % sign changed
        if stat(i) < 0
            sign = 'Negative';
        else
            sign = 'Positive';
        end
    end
    
    fprintf('%s & %s & %s & %d & %.3f & %d %d %d \\\\\n', sign, new_region{i}, '??', extent(i), stat(i), x, y, z);

end

save('shit.mat');
new_region'
