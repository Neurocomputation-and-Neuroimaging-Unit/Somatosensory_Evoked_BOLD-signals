

working_dir = '...\SEBs\ROI_analysis';
addpath('...\spm12');

result_dir = '...\SEBs\Data\2nd_level_GLM_Fcon';

file_3D = niftiread([result_dir filesep 'spmT_0002.nii']);
info = niftiinfo([result_dir filesep 'spmT_0002.nii']);

atlas_path = '...\atlasses\HCPex-main\HCPex_v1.1';
regions = load([atlas_path filesep 'HCPex_2mm_List.mat']).ROI;
atlas_name = 'myMNI_HCPex_2mm.nii';


mask_dir = '...\SEBs\Data\2nd_level_GLM_Fcon\ttest_significant.nii';


mask_3d = niftiread(mask_dir);
atlas = niftiread([atlas_path filesep atlas_name]);


included_regions = {};
max_vox = [];
percent_covered = [];
X_Y_Z = [];
A_B_C = [];
peak_F = [];
peak_p = [];
counter = 0;

%% loop
for r = [27:30 50:52 207:210 230:232] %somatosensory ROIs

    overlap = zeros(79,95,79);
    otherwise_empty = zeros(79,95,79);
    overlap(atlas == r) = 1;
    overlap = overlap.*cast(mask_3d>0, 'double');

    num_vox = sum(overlap>0, 'all');
    this_percent_covered = (num_vox*100)/sum(atlas==r, 'all');

    % for_visuals = cast(num_vox, 'single');


    if num_vox < 1
        continue;
    else

        max_vox = [max_vox; num_vox];
        percent_covered = [percent_covered; this_percent_covered];

        otherwise_empty(overlap>0) = file_3D(overlap>0);

        [peakVal, peakVox] = max(otherwise_empty, [], 'all');
        [X, Y, Z] = ind2sub(size(file_3D), peakVox);
        A_B_C = [A_B_C; X, Y, Z];

        included_regions = {included_regions{:}, regions(r).Nom_L};
        coordsMNI = ([X-1, Y-1, Z-1, 1] * info.Transform.T);
        X_Y_Z = [X_Y_Z; coordsMNI(1:3)]; % spm starts indexing at 0, apparently
        peak_F = [peak_F; peakVal];
        peak_p = [peak_p; 1-fcdf(peakVal,1,2130)];

    end
end

region_names = included_regions';
result_table = table(region_names, max_vox, percent_covered, X_Y_Z, peak_F, peak_p, A_B_C)