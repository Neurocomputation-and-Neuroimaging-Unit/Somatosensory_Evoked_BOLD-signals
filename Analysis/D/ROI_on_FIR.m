% ROI analysis for SEBs


working_dir = '...\SEBs\ROI_analysis';
addpath('...\spm12');

result_dir = 'D:\SEBs\Data\2nd_level_4er_72_bins_1-way-ANOVA';

atlas_path = '...\atlasses\HCPex-main\HCPex_v1.1\';
regions = load([atlas_path filesep 'HCPex_2mm_List.mat']).ROI;
atlas_name = 'myMNI_HCPex_2mm.nii';

reslice = 0;
template = [result_dir filesep 'spmT_0001.nii'];

write_mask = 1;
mask_dir = '...\permutationTesting\3D_masks_s4_positive';


%% reslice
if reslice == 1

    images = {template, [atlas_path filesep atlas_name]};

    matlabbatch{1}.spm.spatial.realign.write.data = images';

    matlabbatch{1}.spm.spatial.realign.write.roptions.which = [1 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 0;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'myMNI_';

    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
    clear matlabbatch

    % update prefix
    atlas_name = ['myMNI_' atlas_name];
end

%% prepare
% create single 3D mask by OR operation on all the separate 3D masks
file_4D = zeros([79, 95, 79, 72]);

if write_mask == 1
    images = {};
    expression = '(';
end

counter = 0;
for b = 73:72*2

    counter = counter + 1;

    if write_mask == 1
        images = {images{:}, [mask_dir filesep 'mask_3D_3clusters_bin-' num2str(counter) '.nii']};
        expression = [expression 'i' num2str(counter) ' '];
        if counter<72
            expression = [expression '+ '];
        end
    end

    if b < 10
        f_str = [result_dir filesep 'spmT_000' num2str(b) '.nii'];
    elseif b < 100
        f_str = [result_dir filesep 'spmT_00' num2str(b) '.nii'];
    else
        f_str = [result_dir filesep 'spmT_0' num2str(b) '.nii'];
    end
    file_4D(:,:,:,counter) = niftiread(f_str);

    m_str = [mask_dir filesep 'mask_3D_3clusters_bin-' num2str(counter) '.nii'];
    mask_4D(:,:,:,counter) = niftiread(m_str);
end
info = niftiinfo(f_str);

if write_mask == 1

    expression = [expression ')>0'];

    matlabbatch{1}.spm.util.imcalc.input = images';
    matlabbatch{1}.spm.util.imcalc.output = 'combined_3Dmask';
    matlabbatch{1}.spm.util.imcalc.outdir = {working_dir};
    matlabbatch{1}.spm.util.imcalc.expression = expression;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
    clear matlabbatch
end

mask_3d = niftiread([working_dir filesep 'combined_3Dmask.nii']);
atlas = niftiread([atlas_path filesep atlas_name]);

included_regions = {};
max_vox = [];
percent_covered = [];
X_Y_Z = [];
A_B_C = [];
peak_time = [];
first_time = [];
last_time = [];
peak_F = [];
peak_p = [];
counter = 0;

%% loop
for r = [27:30 50:52 207:210 230:232] % somatosensory ROIs

    otherwise_empty = zeros(79,95,79,72);
    overlap = zeros(79,95,79);
    overlap(atlas == r) = 1;
    overlap = overlap.*cast(mask_3d>0, 'double');
    overlap4D = repmat(overlap, [1,1,1,72]);
    overlap4D = overlap4D.*cast(mask_4D>0, 'double');

    num_vox_per_bin = squeeze(sum(overlap4D>0, [1,2,3]));
    [this_max_vox, vis_max] = max(num_vox_per_bin);
    this_percent_covered = (this_max_vox*100)/sum(atlas==r, 'all');

    for_visuals = cast(reshape(overlap4D(:, :, :,vis_max)>0, [79,95,79]), 'single');


    if this_max_vox < 1
        continue;
    else

        max_vox = [max_vox; this_max_vox];
        percent_covered = [percent_covered; this_percent_covered];

        otherwise_empty(overlap4D>0) = file_4D(overlap4D>0);

        [peakVal, peakVox] = max(otherwise_empty, [], 'all');
        [X, Y, Z, T] = ind2sub(size(file_4D), peakVox);
        A_B_C = [A_B_C; X, Y, Z];

        [~,~,~,T1] = ind2sub(size(file_4D), find(otherwise_empty,1,'first'));
        [~,~,~,T2] = ind2sub(size(file_4D), find(otherwise_empty,1,'last'));

        included_regions = {included_regions{:}, regions(r).Nom_L};
        coordsMNI = ([X-1, Y-1, Z-1, 1] * info.Transform.T);
        X_Y_Z = [X_Y_Z; coordsMNI(1:3)]; % spm starts indexing at 0, apparently
        first_time = [first_time; (T1-1)*0.18];
        peak_time = [peak_time; (T-1)*0.18];
        last_time = [last_time; (T2-1)*0.18];
        peak_F = [peak_F; peakVal];
        peak_p = [peak_p; 1-fcdf(peakVal,1,2130)];

    end
end

region_names = included_regions';
result_table = table(region_names, max_vox, percent_covered, X_Y_Z, first_time, peak_time, last_time, peak_F, peak_p, A_B_C)

