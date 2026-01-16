function C01_first_level_batch_BIDS_function(only_subject)

clear runs

%%
SPM_path = '.../spm12';
addpath(SPM_path)
addpath(genpath('.../hMRI-toolbox'))

%data source directory
src_dir      = '.../SEBs/Data';
logDir='.../SEBs/Logs';


%subject identifiers if all subjects are to be included
%%%% to do: SJs to sub
cd(src_dir)
clear SJs
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

display('Subjects found:')
SJs % analysis for these subjects
display('Subjects to exclude:')
excludeSJ = [1:33]
disp(only_subject)
excludeSJ(only_subject) = []; % to make parallelization on hpc and usage as standalose easier

ex_sub = only_subject;
excludeRuns = [];

SEBs = 1; % if 1, work on reordered data!


%% selection of analysis steps (1-5) to be performed
analysis_switch = [ 1 2 3 8 ] % 1 2 3 4 5

start_prefix='s4wSEBr'; %eg. s8wra % s8wSEBsite1rb



%session & run identifiers
sessNum = 0;
if exist([src_dir filesep SJs{1} filesep 'ses-1'])==7 
    cd([src_dir filesep SJs{1}])
    sd = dir('ses*');
    sessNum = length(sd);
    for sess = sessNum
        sessions(1, sess) = {sd(sess).name};
    end
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep sessions{1} filesep 'func']);
        rd = dir('sub*.nii');
        for r = 1:length(rd)
            if ismember(r, excludeRuns)
                continue;
            else
                runs(sb, r) = {rd(r).name};
            end
        end
    end
else
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep 'func']);
        rd = dir('sub*.nii');
        counter = 0;
        for r = 1:length(rd)
            if ismember(r, excludeRuns)
                continue;
            else
                counter = counter +1;
                runs(sb, counter) = {rd(r).name};
            end
        end
    end
end

%anatomy identifier
ana=['anat'];
if thisMulti == '1'
    example_runs = runs{ex_sub, multi(ex_sub,:)==1}
else
    example_runs = runs{ex_sub, multi(ex_sub,:)==4}
end

nifti_files = dir(fullfile([src_dir filesep SJs{ex_sub} filesep 'func'], example_runs)); 
anat_files = dir(fullfile([src_dir filesep SJs{ex_sub} filesep ana], ['sub-', '*T1w.nii'])); %look for all anat nifti files

%now we get the data from the json file
json_files = (dir(fullfile([src_dir filesep SJs{ex_sub} filesep 'func'], [example_runs(1:(length(example_runs)-4)), '*json']))); %extract all json files, althoguh they should have the same info
%because for some datasets (flicker and Ganzfeld) we have json files for
%each nift file and these are named differently, we have to check if the
%first command returns an empty structure. If yes, it means the json files
%have a different naming, starting with subject
if isequal(size(json_files), [0, 1])
    json_files = (dir(fullfile([src_dir filesep SJs{ex_sub} filesep 'func'], example_runs)));
end

json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from
TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = length(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order
slice_order = y';

%now get the same info from nifti header
nifti_file_metadata = [nifti_files(1).folder, filesep, nifti_files(1).name];
info = niftiinfo(nifti_file_metadata);
TR_nifti = info.PixelDimensions(4);
n_slices_nifti = info.ImageSize(3);
vox_size=repmat(info.PixelDimensions(1),1,3);

%compare json and nifti header
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti")
end
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti")
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1: extract/set onsets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2:  1st level glm
% must include variable 'onsets' (cell with runs x conditions) including onset-times in sec
shuff = 0;
stim_sites = 1;

condnames = {'stim'};
duration = 0; % epoch duration; for single events set 0

fmri_t = 36; % Microtime resolution; If you have performed slice-timing correction, change this parameter to match the number of slices specified there; otherwise, set default 16
fmri_t0	= 1;
slice_order_stc2 = 1:36;

beta_dir = ['1st_level_detrended_notupsampled']; % folder that will contain the created job.mat file and SPM file

corr_par = 'none'; % default: 'AR(1)', new: 'FAST'

tr = 0.06;

hpf      = 128; % High-pass filter cut-off; default 128 sec --> option to modulate hpf to make shorter
% include multiple regressors (1=yes)
% if 1 (yes), 'hm' and/or 'cc' will be appended to outputfolder_1st
hm=0;   % head motion parameters from realignment (step 4 in B0_preprocessing)
cc=0;   % CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
% do you want to normalize the realigned data too, to use masks in MNI-space for ROI-based decoding?
% if yes (1) nomalization will be initiated before the construction of the glm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 3:  contrasts for GLM
analysisfolder = [beta_dir]; % '_hm' '_cc'
cnames = {'stim'};

cvecs = {[1]};

del=1; % Delete existing contrasts (1=yes)
% were multiple regressors included in 1st level (step 1)?
n_hm=0;   % number of head motion parameters from realignment (step 4 in B0_preprocessing)
n_cc=0;   % number of CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
% if >0, "zeros" will be appended in design matrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 4:  compute FIR

just_contrasts = 0;
SEBs_FIR = 1;
scaling = 'None'; % 'None' or 'Scaling' scales the data
if SEBs_FIR == 1
    prefix_FIR = 's4wSEBr';
    tr_fir = 0.06;
    bin_FIR = 0.18; %in sec --> thats x data points
else
    prefix_FIR = 's8wr';
    tr_fir = tr;
    bin_FIR = tr;
end
total_FIR = 16.02; %in secs needs to represent every TR ideally
condnamesFIR = {'Stim'};
outputfolder_1st_FIR = ['FIR_REO_detrended_notupsampled_' num2str(round(total_FIR/bin_FIR)) '_bins'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% begin
currPrefix = start_prefix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = analysis_switch

    switch n

        %% extract onsets
        case 1

            for s = 1:numel(SJs)

                if ismember(s, excludeSJ)
                    continue;
                else
                    cd(logDir)
                    oruns=dir(['LogFile_subject-' num2str(s) '_*.tsv']);


                    if ~isempty(runs)
                        identifier = 0;
                        for r=1:size(runs,2)


                            this_log = tdfread([logDir filesep oruns(r).name]);

                            % cond 1:
                            if SEBs == 1
                                onsets{s,r,1} = 0;
                            else
                                onsets{s,r,1} = round(this_log.TriggerTime(this_log.CatchTrial == 0)); %%%
                            end

                        end
                    end
                end
            end

            %% glm
        case 2

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    display(['Step 3, 1st level glm: ' SJs{sj} ])
                    subj_dir = fullfile(src_dir, SJs{sj});

                    C12_glm_1stLevel(SJs, sj, subj_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, duration, hm, cc, stim_sites, corr_par);

                end
            end


            %% contrasts
        case 3

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    SJ_dir = [src_dir filesep SJs{sj}];
                    C13_contrast_1stLevel(SJ_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs);

                end

            end

            %% 2nd level
        case 4

            C14_1st_level_FIR(src_dir, SJs, excludeSJ, outputfolder_1st_FIR, prefix_FIR, total_FIR, bin_FIR, condnamesFIR, tr_fir, runs, logDir, SEBs_FIR, just_contrasts, scaling);

    end
end
end

