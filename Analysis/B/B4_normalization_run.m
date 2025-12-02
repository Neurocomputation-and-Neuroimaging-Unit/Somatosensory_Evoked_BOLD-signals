function B7_normalization_run(data_dir, struct_dir, sj, runs, vox_size, currPrefix, r)


['^' currPrefix runs{sj,r}]

spm('defaults','fmri');
spm_jobman('initcfg');

warning off

f1 = spm_select('List', struct_dir, '^y.*\.nii$');
numVols = size(f1,1);
parameter_file = cellstr([repmat([struct_dir filesep], numVols, 1) f1 ]);


f2 = spm_select('List', data_dir, ['^' currPrefix runs{sj,r}])
numVols2 = size(f2,1);

V=spm_vol([repmat([data_dir filesep], numVols2, 1) f2 ',:']);
for i=1:size(V,1)
    functional_imgs{i} = [repmat([data_dir filesep], numVols2, 1) f2 repmat([',' int2str(i)] , numVols2, 1)];
end


%--------------------------------------------------------------------------
% --------------------------- Normalize  ----------------------------
%--------------------------------------------------------------------------
%%
matlabbatch{1}.spm.spatial.normalise.write.subj.def = parameter_file;

%%
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = functional_imgs';

%%
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70;78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox_size;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;


fprintf(['Normalizing ', num2str(sj), '\n']);
spm_jobman('run', matlabbatch)
clear jobs
