function B3_calculate_SEBS(func_folder,func_name,info_name, n_slices, DesignMat, multi_fact, slice_order, wait_vol, column_idices, stim_sites, new_tr)

for s = 1:stim_sites % in case of mult conditions
    %% sort slices into volumes
    % "create" empty volume for first slice
    func_file = spm_select('List',func_folder, func_name);
    info_file = spm_select('List',func_folder, info_name);

    DesignMat(DesignMat(:,6)==1, :) = [];

    func = load_untouch_nii([func_folder filesep func_file]);

    func = func.img;

    volume_0 = load_untouch_nii([func_folder filesep func_file]);

    volume_1 = load_nii_hdr([func_folder filesep info_file]);
    volume_0.hdr = volume_1;
    volume_0.hdr.dime.dim(5)=(wait_vol-2)*(n_slices/multi_fact); %%% to be on the safe side
    volume_0.hdr.dime.pixdim(5)=new_tr; %new TR in s

    vol_0 = func*0;
    vol_0(:,:,:,((wait_vol-1)*(n_slices/multi_fact)+1):end) = [];

    % fill volume with slices that were stim slices
    for v = 1:(wait_vol-2)*(n_slices/multi_fact) % go trough volumes: from stim slice to volumes-2 (to have same amount of data without taking data correlated w stim slices last)

        % find stim_slice: index of acq slice OR slice(s) after
        if stim_sites > 1
            stim_slice_ind = DesignMat((DesignMat(:,3)==s),column_idices)+(v-1);
        else
            stim_slice_ind = DesignMat(:,column_idices)+(v-1);
        end
        retreived_slice_idx = mod(stim_slice_ind,(n_slices/multi_fact)); % this was the xth slice acq in this vol
        retreived_slice_idx(find(~retreived_slice_idx)) = (n_slices/multi_fact);

        % find according volume
        stim_vol_ind = ceil(stim_slice_ind/(n_slices/multi_fact));
        % find stimulated slice within that volume: mapping from slice time
        % to location (dorsal-ventral)
        %%%%%%%%%%%%%% same for multi
        %  put that slice in the volume
        for sl = 1:(n_slices/multi_fact)
            sliceInVol_idx = slice_order(((retreived_slice_idx(sl)-1)*multi_fact + 1):retreived_slice_idx(sl)*multi_fact);
            vol_0(:,:,sliceInVol_idx,v)=func(:,:,sliceInVol_idx,stim_vol_ind(sl));
        end
    end

    volume_0.img = vol_0;
    if stim_sites > 1
        save_untouch_nii(volume_0, [func_folder filesep 'SEBsite' num2str(s) func_name(2:end-4)]);
    else
        save_untouch_nii(volume_0, [func_folder filesep 'SEB' func_name(2:end-4)]);
    end


    clear func volume_0 vol_0 volume_1

end
end
