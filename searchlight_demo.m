r = 1.814;

[mask, Vmask] = load_mask('masks/mask.nii');
Vmask.fname = 'masks/searchlight.nii';

[x, y, z] = ind2sub(size(mask), find(mask)); % binary mask --> voxel indices --> voxel coordinates in AAL2 space

% randomize the order
%
%idx = randperm(length(X));
%x = x(idx);
%y = y(idx);
%z = z(idx);


min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);
min_z = min(z);
max_z = max(z);

newmasks = {};

for i = randi(numel(x))
    
    newmask = zeros(size(mask));
    for newx = floor(x(i) - r) : ceil(x(i) + r)
        if newx < min_x || newx > max_x, continue; end
        for newy = floor(y(i) - r) : ceil(y(i) + r)
            if newy < min_y || newy > max_y, continue; end
            for newz = floor(z(i) - r) : ceil(z(i) + r)
                if newz < min_z || newz > max_z, continue; end
                if ~mask(newx, newy, newz), continue; end
                if (x(i) - newx)^2 + (y(i) - newy)^2 + (z(i) - newz)^2 > r^2, continue; end
                newmask(newx, newy, newz) = 1;
            end
        end
    end
    
    mni = cor2mni([x(i) y(i) z(i)], Vmask.mat);
    fprintf('%d,%d,%d (%d,%d,%d): %d\n', x(i), y(i), z(i), mni, sum(newmask(:)));
    spm_write_vol(Vmask, newmask);
end

check_mask(Vmask.fname);
