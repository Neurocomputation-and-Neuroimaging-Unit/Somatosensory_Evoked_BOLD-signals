function B0_prepro_SEBs_Batch_function(only_subject)

%%
SPM_path = '.../spm12';
addpath(SPM_path)
addpath(genpath('.../hMRI-toolbox'))
addpath(genpath('.../niftiiTools'))


%data source directory
src_dir      = '.../SEBs/Data';

%subject identifiers
cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

excludeRuns = [];

excludeSJ = [1:33]; % zb. unvollständige datensätze: 1 22
%session & run identifiers
excludeSJ(only_subject) = [];

disp(num2str(only_subject))

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
        for r = 1:length(rd)
            if ismember(r, excludeRuns)
                continue;
            else
                runs(sb, r) = {rd(r).name};
            end
        end
    end
end

%anatomy identifier
ana=['anat'];
%unzip
zip_files = dir(fullfile(src_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files)
    for z = 1:size(zip_files, 1)
        gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
    end
end

%% selection of analysis steps (1-12) to be performed
analysis_switch = [1 2 3 4 5 6 7]; % 3 4 7 8
start_prefix= ''; % e.g. s4wSEBdr

%now we get the data from the json file
json_files = (dir(fullfile(src_dir, '**', ['task', '*json']))); %extract all json files, althoguh they should have the same info
%because for some datasets (flicker and Ganzfeld) we have json files for
%each nift file and these are named differently, we have to check if the
%first command returns an empty structure. If yes, it means the json files
%have a different naming, starting with subject
if isequal(size(json_files), [0, 1])
    json_files = (dir(fullfile(src_dir, '**', ['sub-', '*bold.json'])));
end

json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = length(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order

n_slices = n_slices_json;

%% spec steps

% 1 realign

% 2 detrend

% 3 calc SEBs
runs_to_SEB = [1 2 3 4];
log_folder = '.../SEBs/Logs';
column_of_matrix = 4;% which column represents index of slice in log
stim_sites = 1;
new_tr = 0.06; %% carefull: depends on sequence

% 4 coregister
corrPrefix = '';

% 5 segment

% 6 normalize
vox_size=[3 3 3]; % voxel size in mm

% 7 smooth
kernel_size=[4 4 4]; %FWHM kernel size



shuff = 0;
start_at_shuff = 0;


%% execution of jobs

currPrefix=start_prefix;

for n=analysis_switch

    switch n

        case 1 % realign

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    sj_dir = [src_dir filesep SJs{sj}];
                    if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, 1}])
                        display(['Step 1, realignment: ' SJs{sj} ', ' runs{sj, 1}])
                        funcPath = [src_dir filesep SJs{sj} filesep 'func'];
                        run_dir = fullfile(funcPath);
                        counter = 0;
                        for r = 1:size(runs, 2)
                            if ismember(r, excludeRuns)
                                continue;
                            else

                                counter = counter+1;
                                run_files{counter} = spm_select('List',run_dir,['^' currPrefix runs{sj, r}]);
                            end
                        end
                        B1_Realignment_all_runs(sj_dir, run_files);

                    elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                        for ses = 1:sessNum
                            display(['Step 1, realignment: ' SJs{sj} ', ' runs{sj, r}])
                            sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                            run_dir = fullfile(sesPath);
                            for r = 1:size(runs, 2)
                                run_files{r} = spm_select('List',run_dir,['^' currPrefix runs{sj, r}]);
                            end
                            B1_Realignment_all_runs(sesPath, run_files);
                        end
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                    end
                    %end
                end
            end

            currPrefix=['r' currPrefix];
       case 2 % 2 detrend

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    display(['Step 2, Detrending: ' SJs{sj}])
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])

                            funcPath = [src_dir filesep SJs{sj}];
                            run_dir = fullfile(funcPath, 'func');
                            f = spm_select('List',run_dir, ['^' currPrefix runs{sj, r}]);
                            V=spm_vol([run_dir filesep f(1,:)]);
                            files={};
                            for i=1:(size(V,1))
                                files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
                            end
                            fileset{r}=char(files);

                        elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                sessPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                run_dir = fullfile(sessPath, 'func');
                                f = spm_select('List',run_dir, ['^' currPrefix runs{sj, r}]);
                                V=spm_vol([run_dir filesep f(1,:)]);
                                files={};
                                for i=1:(size(V,1))
                                    files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
                                end
                                fileset{r}=char(files);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                    if exist('fileset')
                        B2_detrending_lmgs(fileset);
                        clear fileset files
                    end
                    display([SJs{sj} ', is done'])
                end
            end
            currPrefix=['d' currPrefix];

        case 3 % SEBs

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    for r = runs_to_SEB
                        example_run = runs{sj, r};

                        json_files = (dir(fullfile(src_dir, '**', [example_run(1:(length(example_run)-4)), '*json']))); %extract all json files, althoguh they should have the same info
                        if isequal(size(json_files), [0, 1])
                            json_files = (dir(fullfile(src_dir, '**', example_run)));
                        end
                        json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from
                        TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
                        slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
                        n_slices_json = length(slice_timing); %compute number of slices from slice timing
                        [~,y]= sort(slice_timing); %compute slice order

                        slice_order = y';
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 3, SEB ' SJs{sj} ', ' runs{sj, r}])
                            run_dir = fullfile(src_dir, SJs{sj}, 'func');
                            cd(log_folder)
                            logs=dir(['DesignMat_subject-' num2str(sj) '_run-' num2str(r) '*']);
                            DesignMat = load([log_folder filesep logs(end).name]);
                            DesignMat = DesignMat.DesignMat;
                            multi_fact = str2num(logs(end).name(regexp(logs(end).name, 'multi-')+6));
                            wait_vol = 8*multi_fact; % or as per script

                            B3_calculate_SEBS(run_dir,['^' currPrefix runs{sj, r}],['^' currPrefix runs{sj, r}], n_slices, DesignMat, multi_fact, slice_order, wait_vol, column_of_matrix, stim_sites, new_tr);

                        elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                display(['Step 3, SEB ' SJs{sj} ', ' runs{sj, r}])
                                sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                run_dir = fullfile(sesPath, 'func');
                                cd(log_folder)
                                logs=dir(['DesignMat_subject-' num2str(sj) '_run-' num2str(r) '_slices-' num2str(n_slices) '_multi-' num2str(multi_fact) '*']);
                                DesignMat = load([log_folder filesep logs(end).name]).DesignMat;
                                B3_calculate_SEBS(run_dir,['^' currPrefix runs{sj, r}], n_slices, DesignMat, multi_fact, slice_order, wait_vol, column_of_matrix, stim_sites, new_tr);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end

            currPrefix=['SEB' currPrefix];

            
        case 4 % coregister

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    %                     for r = 1:size(runs, 2)
                    if exist([src_dir filesep SJs{sj} filesep 'func'])
                        display(['Step 4, coregistration: ' SJs{sj}])
                        funcPath = [src_dir filesep SJs{sj}];
                        func_dir        = fullfile(funcPath, 'func');
                        struct_dir      = fullfile(funcPath, ana);

                        B4_coregister_est(func_dir, struct_dir, sj, '^s.*\.nii', runs, corrPrefix);

                    elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func'])
                        for ses = 1:sessNum
                            display(['Step 4, coregistration: ' SJs{sj}])
                            sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                            func_dir        = fullfile(sesPath, 'func');
                            struct_dir      = fullfile(sesPath, ana);
                            B4_coregister_est(func_dir, struct_dir, sj, '^s.*\.nii', runs, corrPrefix);
                        end
                    else
                        display('###########################################################')
                        display(['############### ' SJs{sj} 's functional data do not exsist ###########'])
                    end
                    %                     end
                end
            end

        case 5 % segment

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                elseif exist([src_dir filesep SJs{sj} filesep ana])==7
                    display(['Step 5, segmentation: ' SJs{sj}])
                    struct_dir = fullfile(src_dir, SJs{sj}, ana);
                    B5_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');
                elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep ana])==7
                    for ses = 1:sessNum
                        sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep ana];
                        display(['Step 5, segmentation: ' SJs{sj}])
                        struct_dir = fullfile(sesPath);
                        B5_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');
                    end
                else
                    display('###########################################################')
                    display(['############### ' SJs{sj} ', ' ana ' does not exsist ###########'])
                end
            end


        case 6 % norm

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else

                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func'])
                            display(['Step 6, normalization: ' SJs{sj}])
                            funcPath = [src_dir filesep SJs{sj}];
                            struct_dir = fullfile(funcPath, ana);
                            data_dir = fullfile(funcPath, 'func');

                            if shuff > 0 && regexp(currPrefix, 'SEB')
                                for sh = start_at_shuff:(start_at_shuff+shuff-1)
                                    B6_normalization_run(data_dir, struct_dir, sj, runs, vox_size, ['SEBshuffle_' num2str(sh) currPrefix(regexp(currPrefix, 'SEB')+3:end)], r);
                                end
                            end

                            B6_normalization_run(data_dir, struct_dir, sj, runs, vox_size, currPrefix, r);

                        elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func'])
                            for ses = 1:sessNum
                                display(['Step 6, normalization: ' SJs{sj}])
                                sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                func_dir        = fullfile(sesPath, 'func');
                                struct_dir      = fullfile(sesPath, ana);
                                B6_normalization_run(data_dir, struct_dir, sj, runs, vox_size, currPrefix);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end

            currPrefix=['w' currPrefix];

        case 7 %  smooth

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 7, smoothing: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj}];
                            run_dir = fullfile(funcPath, 'func');

                            if shuff > 0 && regexp(currPrefix, 'SEB')
                                for sh = start_at_shuff:(start_at_shuff+shuff-1)
                                    runs_shuff = string(strjoin(['SEBshuffle_' num2str(sh) currPrefix(regexp(currPrefix, 'SEB')+3:end) runs(sj, r)], ''));
                                    B7_smoothing_run(run_dir, SJs{sj}, strjoin(['^' currPrefix(1:regexp(currPrefix, 'SEB')-1) runs_shuff], ''),kernel_size);
                                end
                            end
                            B7_smoothing_run(run_dir, SJs{sj}, ['^' currPrefix runs{sj, r}],kernel_size);
                            display([SJs{sj} ', ' runs{sj,r} ' is done'])
                        elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                display(['Step 7, smoothing: ' SJs{sj} ', ' runs{sj, r}])
                                sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                                run_dir = fullfile(sesPath);
                                B7_smoothing_run(run_dir, SJs{sj}, ['^' currPrefix runs{sj, r}],kernel_size);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end

            currPrefix=['s' num2str(unique(kernel_size)) currPrefix];


    end

end

end
