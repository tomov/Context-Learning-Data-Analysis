function [regionname, regionidx] = bspm_getregionnames(xyz, atlaslabels, atlas0, XYZmm0)
    if size(xyz,1)~=3, xyz = xyz'; end
    regionname  = repmat({'Location not in atlas'}, size(xyz, 2), 1);
    regionidx   = zeros(size(xyz, 2), 1);
    for i = 1:size(xyz,2)
        regionidx(i) = atlas0(spm_XYZreg('FindXYZ', xyz(:,i), XYZmm0));
        if regionidx(i)
            regionname{i} = atlaslabels.label{atlaslabels.id==regionidx(i)};
        end
    end
