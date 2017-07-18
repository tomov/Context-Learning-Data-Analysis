function Searchlight = rdms_get_searchlight(data, metadata, which_rows, x, y, z, r)

% Compute the searchlight neural RDMs
%
% INPUT:
% data, metadata = subject data and metadata as output by load_data
% which_rows = which rows (trials) to include
% x, y, z = vectors of voxel coordinates in group-level mask coordinate space
% r = sphere radius
%
% OUTPUT:
% Searchlight = struct array of RDMs

searchlight_idx = 0;

disp('Computing searchlight RDMs...');
tic

[mask, Vmask] = load_mask('masks/mask.nii');
Vmask.fname = 'masks/searchlight.nii'; % !!!! IMPORTANT! in case we save it

[all_x, all_y, all_z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space
min_x = min(all_x);
max_x = max(all_x);
min_y = min(all_y);
max_y = max(all_y);
min_z = min(all_z);
max_z = max(all_z);


for event = {'trial_onset', 'feedback_onset'}
    event = event{1};

    % Load neural data for given event
    %
    whole_brain_betas = get_betas('masks/mask.nii', event, data, metadata);

    % For each voxel, build a mask that's a sphere around it
    % Make sure to only include voxels from the brian
    %
    for i=1:length(x)

        % Build spherical mask
        %
        sphere_mask = zeros(size(mask));
        for newx = x(i) - r : x(i) + r
            if newx < min_x || newx > max_x, continue; end
            for newy = y(i) - r : y(i) + r
                if newy < min_y || newy > max_y, continue; end
                for newz = z(i) - r : z(i) + r
                    if newz < min_z || newz > max_z, continue; end
                    if ~mask(newx, newy, newz), continue; end
                    if (x(i) - newx)^2 + (y(i) - newy)^2 + (z(i) - newz)^2 > r^2, continue; end
                    sphere_mask(newx, newy, newz) = 1;
                end
            end
        end

        % Get sphere center coordinates in MNI space
        %
        mni = cor2mni([x(i) y(i) z(i)], Vmask.mat);

        sphere_betas = get_betas_submask(sphere_mask, whole_brain_betas);
        [sphereRDMs, avgSphereRDM] = compute_rdms(sphere_betas, 'cosine', data, metadata, which_rows);
        searchlight_idx = searchlight_idx + 1;
        Searchlight(searchlight_idx).RDMs = sphereRDMs;
        Searchlight(searchlight_idx).RDM = avgSphereRDM;
        Searchlight(searchlight_idx).name = ['sphere_', sprintf('%d_%d_%d', mni), '_', event(1)];
        Searchlight(searchlight_idx).color = [0 1 0];
        Searchlight(searchlight_idx).center = [x(i) y(i) z(i)];
        Searchlight(searchlight_idx).radius = r;
        Searchlight(searchlight_idx).center_mni = mni;
    end

end

disp('Computed searchlight RDMs.');
toc
