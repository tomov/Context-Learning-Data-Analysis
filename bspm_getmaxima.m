function PEAK                    = getmaxima(Z, XYZ, M, Dis, Num)
[N,Z,XYZ,A,L]       = spm_max(Z,XYZ);
XYZmm               = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
npeak   = 0;
PEAK    = [];
while numel(find(isfinite(Z)))
    %-Find largest remaining local maximum
    %------------------------------------------------------------------
    [U,i]   = max(Z);            %-largest maxima
    j       = find(A == A(i));   %-maxima in cluster
    npeak   = npeak + 1;         %-number of peaks
    extent  = N(i);
    PEAK(npeak,:) = [extent U XYZmm(:,i)']; %-assign peak
    %-Print Num secondary maxima (> Dis mm apart)
    %------------------------------------------------------------------
    [l,q] = sort(-Z(j));                              % sort on Z value
    D     = i;
    for i = 1:length(q)
        d = j(q(i));
        if min(sqrt(sum((XYZmm(:,D)-repmat(XYZmm(:,d),1,size(D,2))).^2)))>Dis
            if length(D) < Num
                D          = [D d];
                npeak   = npeak + 1;         %-number of peaks
                PEAK(npeak,:) = [extent Z(d) XYZmm(:,d)']; %-assign peak
            end
        end
    end
    Z(j) = NaN;     % Set local maxima to NaN
end
