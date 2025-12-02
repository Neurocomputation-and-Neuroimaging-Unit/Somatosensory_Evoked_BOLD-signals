function C31_1st_level_FIR_shuff(src_dir, SJs, excludeSJ, fir_out, shuffles_2nd, total_FIR, bin_FIR)


% seed rng as randomly as possible
rng(shuffles_2nd(end)*randi(100));


for sh_2nd = shuffles_2nd

    %% do 2nd level here
    % target directory that will contain the created job.mat file

    tgt_dir_2nd      = [src_dir filesep '2nd_levels_s4_permuted_noDep_real' filesep '2nd_level_bins-' num2str(round(total_FIR/bin_FIR)) '_shuffle1-' num2str(sh) '_shuffle2-' num2str(sh_2nd) '_wANOVA'];

    % Create tgt_dir
    if ~exist(tgt_dir_2nd, 'dir')
        mkdir(tgt_dir_2nd)
    end

    matlabbatch{1}.spm.stats.factorial_design.dir = {tgt_dir_2nd};

    % create SPM style file list for model specification.
    counti = 0;
    for s = 1:numel(SJs)

        if ismember(s, excludeSJ)
            continue;
        else

            counti = counti + 1;
            fNames = [];
            for con = randperm(round(total_FIR/bin_FIR)) %%%%%
                if con < 10
                    filt            = [['con_000' num2str(con)] '*\.nii$'];
                else
                    filt            = [['con_00'  num2str(con)] '*\.nii$'];
                end

                temp_dir      = fullfile(src_dir, SJs{s}, [fir_out]); 
                cd(temp_dir)
                f             = spm_select('List', temp_dir, filt);
                fs            = cellstr([repmat([temp_dir filesep], 1, 1) f repmat(',1', 1, 1)]);
                fNames          = [fNames; fs];
            end
            matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(counti).scans = fNames;
            matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(counti).conds = (1:round(total_FIR/bin_FIR))';
        end
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
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(tgt_dir_2nd, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch



end