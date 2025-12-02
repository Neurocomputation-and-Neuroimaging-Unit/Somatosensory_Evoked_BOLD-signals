function C21_2ndLevel_OneSampleTest(src_dir, SJin, outputfolder, dir_1st, cnames)

% This function performs a second level statistical analysis by performing
% one-sample t/F-tests over .con images

% target directory that will contain the created job.mat file
tgt_dir      = [src_dir filesep outputfolder]

% Create tgt_dir
if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end

for dir = 1:length(dir_1st)
    % cycle over contrasts
    for con = 1:numel(cnames)
        % obtain the single subject contrast image filenames
        fNames = [];
        for sj = 1:numel(SJin)
            % create SPM style file list for model specification.
            if con < 10
                filt            = [['con_000' num2str(con)] '*\.nii$'];
            else
                filt            = [['con_00'  num2str(con)] '*\.nii$'];
            end

            temp_dir      = fullfile(src_dir, SJin{sj}, dir_1st{dir})
            cd(temp_dir)
            f             = spm_select('List', temp_dir, filt);
            fs            = cellstr([repmat([temp_dir filesep], 1, 1) f repmat(',1', 1, 1)]);
            fNames          = [fNames; fs];

        end
    end

    cur_dir = tgt_dir;
    if ~exist(cur_dir, 'dir')
        mkdir(cur_dir)
    end

    % create one sample t-test GLM
    % -------------------------------------------------------------------------
    matlabbatch{1}.spm.stats.factorial_design.dir                       = {cur_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = cellstr(fNames);
    matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em                = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;

    %%  Model Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(cur_dir, 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    %% T-Contrast Specification
    % one-sample t-test contrast specification

    connames = 'F-Contrast';
    convecs  = 1;

    % Number of contrasts
    numCons  = 1;

    % Allocate SPM.mat file
    %
    matlabbatch{3}.spm.stats.con.spmmat = {fullfile(cur_dir, 'SPM.mat')};
    % Cycle over contrast specifications
    for c = 1:numCons
        % Allocate F-contrast structure
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name    = connames;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = convecs;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
    end
    % Delete existing contrasts (1=yes)
    matlabbatch{3}.spm.stats.con.delete = 1;

    %% Run the job ['Con' num2str(con) '_' cnames{con}]
    fprintf(['Computing 2nd Level Contrast ' num2str(con) ': ' cnames{con} '\n'])
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);

    % clear job variable
    clear matlabbatch
end
end
