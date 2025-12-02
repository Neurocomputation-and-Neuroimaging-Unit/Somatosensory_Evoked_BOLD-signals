function C22_2ndLevel_1wayANOVA(src_dir, SJin, outputfolder, dir_1st, con_images)

%% 
% This function performs a second level statistical analysis 

% target directory that will contain the created job.mat file
tgt_dir      = [src_dir filesep outputfolder];

% Create tgt_dir
if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end

matlabbatch{1}.spm.stats.factorial_design.dir = {tgt_dir};

for dir = 1:length(dir_1st)
        % create SPM style file list for model specification.

        for sj = 1:numel(SJin)

            fNames = [];

            for con = con_images
                if con < 10
                    filt            = [['con_000' num2str(con)] '*\.nii$'];
                else
                    filt            = [['con_00'  num2str(con)] '*\.nii$'];
                end

                temp_dir      = fullfile(src_dir, SJin{sj}, dir_1st{dir});
                cd(temp_dir)
                f             = spm_select('List', temp_dir, filt);
                fs            = cellstr([repmat([temp_dir filesep], 1, 1) f repmat(',1', 1, 1)]);
                fNames          = [fNames; fs];
            end

            matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(sj).scans = fNames;
            matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(sj).conds = con_images';

        end
     
end


matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 0; 
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 0; 
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


%%  Model Estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(tgt_dir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch

end