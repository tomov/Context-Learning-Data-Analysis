%------------ FreeSurfer -----------------------------%
fshome = getenv('FREESURFER_HOME');
fsmatlab = sprintf('%s/matlab',fshome);
if (exist(fsmatlab) == 7)
    path(path,fsmatlab);
end
clear fshome fsmatlab;
%-----------------------------------------------------%

%------------ FreeSurfer FAST ------------------------%
fsfasthome = getenv('FSFAST_HOME');
fsfasttoolbox = sprintf('%s/toolbox',fsfasthome);
if (exist(fsfasttoolbox) == 7)
    path(path,fsfasttoolbox);
end
clear fsfasthome fsfasttoolbox;
%-----------------------------------------------------%

%------------ SPM 

addpath(genpath(getenv('_HVD_SPM_DIR')));

% do NOT add(genpath(...) -- it makes MATLAB loading way too slow
addpath('/ncf/gershman/Lab/scripts/spm12');
addpath('/ncf/gershman/Lab/scripts/ccnl-fmri');
addpath('/ncf/gershman/Lab/scripts/mfit');
addpath('/ncf/gershman/Lab/scripts/bspmview');
addpath('/ncf/gershman/Lab/scripts/helper');
addpath('/ncf/gershman/Lab/scripts/matlab/ConLearn/libs/rsatoolbox/Engines');
addpath('/ncf/gershman/Lab/scripts/matlab/ConLearn/libs/rsatoolbox/Modules');
%addpath(genpath('/ncf/gershman/Lab'));

