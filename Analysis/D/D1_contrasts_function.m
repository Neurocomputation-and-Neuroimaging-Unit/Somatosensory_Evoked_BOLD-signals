function D1_contrasts_function(this_dir, bins)
% contrasts

matlabbatch{1}.spm.stats.con.spmmat = {[this_dir filesep 'SPM.mat']};

% hardcoded sj number and contrast structure - for now
% change for later versions!
cons = [-eye(bins), -ones(bins,31)/31; eye(bins), ones(bins,31)/31]; 


for b = 1:bins*2
    if b < bins+1
        sign = -1;
    else
        sign = 1;
    end
    bin_nr = mod(b,bins);
    if bin_nr == 0
        bin_nr = bins
    end
    matlabbatch{1}.spm.stats.con.consess{b}.tcon.name = [num2str(bin_nr*sign)];
    matlabbatch{1}.spm.stats.con.consess{b}.tcon.weights = cons(b,:);
    matlabbatch{1}.spm.stats.con.consess{b}.tcon.sessrep = 'none';
end
matlabbatch{1}.spm.stats.con.delete = 1;

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

end

