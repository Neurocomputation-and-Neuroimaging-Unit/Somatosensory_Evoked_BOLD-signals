% the automatic region retreiver using the HCPex atlas
% input: x y z in MNI, top N results
% output: top N regions + according probability

function [top_probN] = get_region_from_anatomy_HCPex(xyzmm)

    atlas_path = '...\atlasses\HCPex-main\HCPex_v1.1';
    addpath(atlas_path)
    regions = load([atlas_path filesep 'HCPex_2mm_List.mat']).ROI;
    atlas = niftiread([atlas_path filesep 'myMNI_HCPex_2mm.nii']);
    info = niftiinfo([atlas_path filesep 'myMNI_HCPex_2mm.nii']);

    abc = [xyzmm', 1] * inv(info.Transform.T);
    A = abc(1)+1;
    B = abc(2)+1;
    C = abc(3)+1;

    int_label = atlas(A,B,C);

    if int_label > 0
        str_label = regions(int_label).Nom_L;
        top_probN = str_label;
    else
        top_probN = 'unassigned';
    end

end
