function C03_2nd_level_shuffle_function(shuffles_2nd)
% to parallelise as many permutations as possible

clear runs

%%
SPM_path = '.../spm12';
addpath(SPM_path)

%data source directory
src_dir      = '.../SEBs/Data';
logDir='.../SEBs/Logs';

cd(src_dir)
clear SJs
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

display('Subjects found:')
SJs % analysis for these subjects
display('Subjects to exclude:')
excludeSJ = [1 22]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# settings:  1st level info

bin_FIR = 0.18; %in sec --> thats x data points

total_FIR = 12.96; %in secs needs to represent every TR ideally
outputfolder_1st_FIR = ['FIR_REO_-s4_' num2str(round(total_FIR/bin_FIR)) '_bins'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C31_1st_level_FIR_shuff(src_dir, SJs, excludeSJ, outputfolder_1st_FIR, shuffles_2nd, total_FIR, bin_FIR);


end
