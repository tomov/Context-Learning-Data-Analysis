function [out, outmat] = bspm_reslice_image(in, ref, int)
    % Most of the code is adapted from rest_Reslice in REST toolbox:
    % Written by YAN Chao-Gan 090302 for DPARSF. Referenced from spm_reslice.
    % State Key Laboratory of Cognitive Neuroscience and Learning
    % Beijing Normal University, China, 100875
    % int:        interpolation, 0=Nearest Neighbor, 1=Trilinear(default)
    if nargin<3, int = 1; end
    if nargin<2, display('USAGE: [out, outmat] = reslice_image(infile, ref, SourceHead, int)'); return; end
    if iscell(ref), ref = char(ref); end
    if iscell(in), in = char(in); end
    % read in reference image
    RefHead = spm_vol(ref);
    mat=RefHead.mat;
    dim=RefHead.dim;
    SourceHead = spm_vol(in);
    [x1,x2,x3] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
    d       = [int*[1 1 1]' [1 1 0]'];
    C       = spm_bsplinc(SourceHead, d);
    v       = zeros(dim);
    M       = inv(SourceHead.mat)*mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
    y1      = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
    y2      = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
    y3      = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
    out     = spm_bsplins(C, y1,y2,y3, d);
    tiny = 5e-2; % From spm_vol_utils.c
    Mask = true(size(y1));
    Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
    Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
    Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));
    out(~Mask) = 0;
    outmat = mat;
