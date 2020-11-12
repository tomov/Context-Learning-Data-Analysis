% compare GLMs using BMS

clear all;

EXPT = context_expt();

r = 10; % mm

masks = isl_create_masks(false, r);

masks = {'masks/IFG_Tri_R.nii'};

%glms = [1 178 179 180 181  187 188 189 190];
glms = [1 178 179 180 181];
%glms = [178 179 180 181];
%glms = [1 178 191 192 193];
%glms = [1 178  179 180 181 191 192 193];
%glms = [1 194:197];

%masks = isl_create_masks_KL(true, r);
%glms = [182 183 184 185  ];

for c = 1:length(masks)
    mask = masks{c};

    mask

    lmes = [];
    bics{c} = [];
    for i = 1:length(glms)
        glmodel = glms(i);
        bic = ccnl_bic(EXPT, glmodel, mask, getGoodSubjects());
        bics{c} = [bics{c} bic];
    end

    lme = -0.5 * bics{c};
    [alpha, exp_r, xp, pxp, bor] = bms(lme);

    pxp

    pxps(c,:) = pxp;
    xps(c,:) = xp;
    bors(c,:) = bor;
    exp_rs(c,:) = exp_r;
end

T = table(masks', pxps, bors)


filename = sprintf('isl_glm_bms.mat', r);
%save(fullfile('mat', filename));

