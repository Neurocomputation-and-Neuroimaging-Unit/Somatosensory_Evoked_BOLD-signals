

clear runs

%%
SPM_path = '.../spm12';
addpath(SPM_path)
addpath(genpath('.../hMRI-toolbox'))

%data source directory
src_dir      = '.../SEBs/Data';

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
excludeSJ = [1 22] % incomplete data sets


%% selection of analysis steps (1-5) to be performed
analysis_switch = [ 1 2 ] % 1 2 3 4 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1:  2nd level after GLM: 1 sample t or F test
con_images=1; % bei 2xN: 1:9,b 2xA: 1:10
outputfolder_2nd = '2nd_level_simple-Fcon';
dir_1st   =  {'1st_level_detrended_notupsampled'}; 
cnames_2nd = {'stim'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2:  2nd level ANOVA 1way after FIR
con_images_A = 1:72; % 18	36	45	54
outputfolder_2nd_A = ['2nd_level_' num2str(length(con_images_A)) '_bins_1-way-ANOVA'];
dir_1st_A   =  {['FIR_REO_detrended_notupsampled_' num2str(length(con_images_A)) '_bins']}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
currPrefix = start_prefix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = analysis_switch

    switch n

        %% 2nd level after GLM
        case 1

            SJin=SJs;
            SJin(excludeSJ)=[];
            C21_2ndLevel_OneSampleTest(src_dir, SJin, outputfolder_2nd, dir_1st, cnames_2nd);

        %% 2nd level after FIR
        case 2
            
            SJin=SJs;
            SJin(excludeSJ)=[];
            C22_2ndLevel_1wayANOVA(src_dir, SJin, outputfolder_2nd_A, dir_1st_A, con_images_A);

    end
end

