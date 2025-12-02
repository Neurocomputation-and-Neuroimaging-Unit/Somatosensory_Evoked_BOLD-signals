function [maxArea, keepers, L, file_4D, CC, info] = D2_cluster_detection_4D(dir_2nd, cluster_thresh_p, n2, min_size, num_bins, this_sign)

% cluster_detection_4D searches for connected components throughout many 3D
% NIFTIs that belong to statistical maps created through F-contrasts

% OUTPUT:
%   maxArea:    Voxel/Time-point count of the biggest cluster in this analysis
%   keepers:    4D map of all the clusters resulting from this analysis
%               under the assigned parameters (see input)
%   L:          4D labels of all clusters resulting from the analysis

% INPUT:
%   dir_2nd:            directory containg all the 2nd level F-maps
%   cluster_thresh_p:   cluster-forming threshold as a p-value
%   n2:                 second degrees of freedom of the F-test
%   min_size:           cluster size cutoff for the output "keepers"
%   num_bins:           number of temporal bins or number of F-maps to
%                       include

% cluster_thresh_p = 0.001; % like uncorrected
n1 = 1; % df 1: each F-map separatly  
% n2 = 2130; % df 2: (31-1)*(72-1) in case of F-contrast
% n2 = 2490 if s = 31 and b = 84
% cluster_thresh_F = tinv(1-cluster_thresh_p,n2);
cluster_thresh_T = tinv(1-cluster_thresh_p,n2); % t-contrasts: one - and one +


file_4D = zeros([79, 95, 79, num_bins]);

if this_sign == -1
    con_nums = 1:num_bins;
else
    con_nums = (num_bins+1):num_bins*2;
end

counter = 0;

for sz = con_nums

    counter = counter + 1;

    if sz < 10
        t_str = [dir_2nd filesep 'spmT_000' num2str(sz) '.nii'];
    elseif sz < 100
        t_str = [dir_2nd filesep 'spmT_00' num2str(sz) '.nii'];
    else
        t_str = [dir_2nd filesep 'spmT_0' num2str(sz) '.nii'];
    end
    file_4D(:,:,:,counter) = niftiread(t_str);

end
info = niftiinfo(t_str);

binary_4D = zeros(size(file_4D));

binary_4D(file_4D>cluster_thresh_T) = 1;

CC = bwconncomp(binary_4D);
L = labelmatrix(CC);

prop = regionprops(CC,"Area");
[maxArea,maxIdx] = max([prop.Area]);

[cluster_sizes, sortingIdx] = sort([prop.Area], 'descend'); % so that B = A(I)
kickIdx = sortingIdx(cluster_sizes<min_size);

idx = setdiff(1:CC.NumObjects,kickIdx);
keepers = cc2bw(CC,ObjectsToKeep=idx);



end