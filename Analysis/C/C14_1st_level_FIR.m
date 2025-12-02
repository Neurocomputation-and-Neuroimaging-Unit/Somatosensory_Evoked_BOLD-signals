function C14_1st_level_FIR(src_dir, SJs, excludeSJ, fir_out, currPrefix, total_FIR, bin_FIR, condnames, TR, runs, log_folder, SEBs, just_contrasts, scaling)

% enter data directory here
data_dir  = src_dir;

%% Loop over Subjects
for sj = 1:numel(SJs)

    if ismember(sj, excludeSJ)
        continue;
    else
        % Composition of the Subject directory
        sj_dir = fullfile(data_dir, SJs{sj});

        % SPM defaults
        spm('defaults','fmri');
        spm_jobman('initcfg');
        % OUTPUT Directory (as subdirectory of the SJ directory)
        tgt_dir = fullfile(sj_dir, fir_out);
        if ~exist(tgt_dir, 'dir')
            mkdir(tgt_dir)
        end

        if just_contrasts == 0

            %********************************************************
            %               Setting FIR parameter
            % *******************************************************
            % Output directory
            jobs{1}.stats{1}.fmri_spec.dir = cellstr(tgt_dir);
            % timing parameters
            jobs{1}.stats{1}.fmri_spec.timing.units     = 'secs';
            jobs{1}.stats{1}.fmri_spec.timing.RT        = TR;
            jobs{1}.stats{1}.fmri_spec.timing.fmri_t    = 36;
            jobs{1}.stats{1}.fmri_spec.timing.fmri_t0   = 1;
            % FIR-SPECIFICATION
            jobs{1}.stats{1}.fmri_spec.bases.fir.length = total_FIR; % some seconds
            jobs{1}.stats{1}.fmri_spec.bases.fir.order  = round(total_FIR/bin_FIR);  % many time-bins
            % Other Specifications
            jobs{1}.stats{1}.fmri_spec.fact              = struct('name', {}, 'levels', {});
            jobs{1}.stats{1}.fmri_spec.volt              = 1;
            jobs{1}.stats{1}.fmri_spec.global            = scaling;
            jobs{1}.stats{1}.fmri_spec.mask              = {''};
            jobs{1}.stats{1}.fmri_spec.cvi               = 'None';

            % ********************************************************
            %            Specify the Desing/Conditions/Onsets
            % ********************************************************

            % get all runIDs from file names in 'func' folder
            nifti_dir = fullfile(sj_dir, 'func');

            % Loop over Sessions, as the logfiles (with onset data) are saved
            % separately per session
            counter = 0;

            for r = 1:size(runs,2)

                counter = counter + 1;
                if SEBs == 1
                    onsets = 0;
                else
                    cd(log_folder);
                    this_logs = dir(['LogFile_subject-' num2str(sj) '_run-' num2str(r) '*.tsv']);
                    this_log = tdfread(this_logs(end).name);
                    onsets = round(this_log.TriggerTime(this_log.CatchTrial == 0));
                end

                % high-pass cut-off
                jobs{1}.stats{1}.fmri_spec.sess(counter).hpf     = 128;
                % Allocation of Data (EPIs/Images) for the current Session
                filename = fullfile(nifti_dir, [currPrefix runs{sj,r}]);
                % use 'expand' to read 4d nifti file
                f = spm_select('expand', filename);
                jobs{1}.stats{1}.fmri_spec.sess(counter).scans   = cellstr(f);

                jobs{1}.stats{1}.fmri_spec.sess(counter).cond(1).name     = condnames{1};
                jobs{1}.stats{1}.fmri_spec.sess(counter).cond(1).onset    = onsets;
                jobs{1}.stats{1}.fmri_spec.sess(counter).cond(1).duration = 0;
                jobs{1}.stats{1}.fmri_spec.sess(counter).cond(1).tmod     = 0;
                jobs{1}.stats{1}.fmri_spec.sess(counter).cond(1).pmod     = struct('name', {}, 'param', {}, 'poly', {});
            end

        end

        % Create model
        fprintf(['Creating GLM\n'])
        spm_jobman('run', jobs);
        clear jobs

        %  Model Estimation
        load(fullfile(tgt_dir, 'SPM.mat'));
        fprintf(['Estimating GLM \n']);
        cd(tgt_dir);
        SPM = spm_spm(SPM);
        clear SPM;

    end

    %% now contrasts
    cnames = condnames;
    cvecs = zeros(round(total_FIR/bin_FIR),1);

    % Allocate SPM.mat file
    matlabbatch{1}.spm.stats.con.spmmat = {[tgt_dir filesep 'SPM.mat']};

    % Number of Runs
    % rundir=dir([subj_dir filesep 'func']);
    numRuns=size(runs,2)/2;       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number of contrasts
    numCons  = numel(cnames);

    % Cycle over contrast specifications
    for c = 1:numCons
        for b = 1:round(total_FIR/bin_FIR)

            convec=cvecs;
            convec(b) = 1;

            % Allocate t-contrast structure
            matlabbatch{1}.spm.stats.con.consess{round(total_FIR/bin_FIR)*(c-1)+b}.tcon.name    = [cnames{c} '_' num2str(b)];
            matlabbatch{1}.spm.stats.con.consess{round(total_FIR/bin_FIR)*(c-1)+b}.tcon.weights  = repmat(convec', 1, numRuns);
            matlabbatch{1}.spm.stats.con.consess{round(total_FIR/bin_FIR)*(c-1)+b}.tcon.sessrep = 'none';
        end
    end

    % Delete existing contrasts (1=yes)
    matlabbatch{1}.spm.stats.con.delete = 1;

    % Run the job
    fprintf(['Computing 1st Level Contrasts\n'])
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);

    % clear job variable
    clear matlabbatch
end

end