% compare GLMs using BMS

clear all;

EXPT = context_expt();

r = 4; % mm

masks = isl_create_masks(false, r);
glms = [178 179 180 181];

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
end

T = table(masks', pxps, xps)


save('isl_glm_bms.mat');

