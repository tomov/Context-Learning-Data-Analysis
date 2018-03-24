% Create a statistical map based on the pre-computed classifier searchlight results by classify_searchlight.m in the rdms/ directory
%

clear all;
close all;

dirname = 'classify/patternnet_1classifier';



%% get model names
%
[data, metadata] = load_data('data/fmri.csv', true, getGoodSubjects());
which_rows = data.which_rows & data.isTrain; % Look at training trials only

%% initialize an empty map
%
[~, V, amap] = load_mask(fullfile('masks', 'spmT_0001.nii'));
amap(:) = NaN; % clear
V.fname = fullfile(dirname, ['searchlight_amap.nii']); % change immediately!

%% load all the searchlight results
%
files = dir(dirname);
for i = 1:length(files)
    file = files(i).name;
    if startsWith(file, 'searchlight_') && endsWith(file, '.mat')
        disp(['Loading ', file]);
        table_T = [];
        load(fullfile(dirname, file));
        
        for j = 1:length(x)
            amap(x(j), y(j), z(j)) = Searchlight(j).accuracy;
        end
    end
end

spm_write_vol(V, amap);

%% visualize
%
struc = fullfile('masks','mean.nii');

bspmview(V.fname, struc);
