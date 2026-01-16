% permutation statistics for 4D data following 1st level FIR and then 2nd
% level ANOVA to correct for multiple comparisons introduced by temporal dimension
% finally: report results in 4D table

% roadmap:
% 1.    calculate contrasts per permutation if necessary
% 2.    calculate distribution for comparison by looping over permuted 2nd
%       levels: keep largest (4D) cluster size --> calculate cut-off
% 3.    calculate 4D clusters for analysis of interest --> apply cut-off
% 4.    extract further info from analysis of interst and create
%       statistical output for all significant clusters:
%           - size
%           - cluster p-value (based on comparison distribution)
%           - peak voxel (x y z)
%           - peak time point (bin / s)
%           - peak F-value
%           - peak z-value


%% input: general
addpath('...'); % working dir
addpath('...\spm12');

aoi_dir = '...\SEBs\Data\2nd_level_4er_89_bins_1-way-ANOVA'; % directory of analysis of interest
perm_dir = '...\SEBs\Data\2nd_levels_detrended_notupsampled_permuted'; % directory of permutation analyses

cluster_forming_p = 0.001;
df2 = 2640; % (sj - 1) * (b - 1) = 30*88=2640
num_bins = 72;
all_clusters = 0; % of all clusters = 0 it only takes the largest, if all clusters > 0, the value is the minimal cluster size to contribute
per_cluster = 0; % number of clusters you want to include, e.g. 6

this_sign = 1; % look at positive (1) or negative (-1) clusters?

%% 1: contrasts for all permutations
only_perm = ['']; % run contrasts only for this permutation

%% 2: creation of distribution from permutations
perm_p_val = 0.05; % top x per-cent
load_dist = ['']; % e.g 'distribution_comparison_100perms.mat'
exp_plot = 0;
plot_result = 1;
redistribute_dist = 0; % 0 = no, 1 = yes, flatten, 2 = yes, take maximum statistics

%% 3: results
within_cluster_stats = 0;
min_k_3D = [10];
save_3D_masks = 0;
% mask uncorrected results with these to have 4D error correction:
mask_save_dir = '...\permutationTesting\3D_masks_s4_positive';

%% put steps here:
analysis_switch = [1 2 3 4]; % e.g. 1 2 3

%% compute:
for a = analysis_switch
    switch a
        %% create F-maps
        case 1

            cd(perm_dir);
            perms = dir('2nd*');
            perm_num = size(perms,1);
            disp([num2str(perm_num) ' permuted 2nd levels found. Calculating Fmaps...'])
            for p = 1:perm_num
                this_dir = [perms(p,1).folder filesep perms(p,1).name];
                if ~isempty(only_perm) && ~strcmp(this_dir, [perms(p,1).folder filesep only_perm])
                    continue;
                else
                    disp(this_dir)
                    try
                        D1_contrasts_function(this_dir, num_bins)
                    catch
                        warning(['Computation of permutation ' num2str(p) ' in folder ' this_dir ' failed!']) % check: was download from hpc complete?
                    end
                end
            end

            %% create distribution
        case 2

            cd(perm_dir);
            perms = dir('2nd*');
            perm_num = size(perms,1);
            disp([num2str(perm_num) ' permuted 2nd levels found. Calculating distribution...'])

            if all_clusters == 0
                size_distribution = zeros(size(perms));
            else
                size_distribution = [];
            end

            for p = 1:perm_num

                this_dir = [perms(p,1).folder filesep perms(p,1).name];
                if mod(p,10) == 0
                    disp(['Currently reading the ' num2str(p) 'th permutation...'])
                end
                if all_clusters == 0
                    [maxArea_p, ~, ~, ~, ~, ~] = D2_cluster_detection_4D(this_dir, cluster_forming_p, df2, 1, num_bins, this_sign);
                    size_distribution(p,1) = maxArea_p;
                elseif per_cluster == 0
                    [~, ~, ~, ~, CC, ~] = D2_cluster_detection_4D(this_dir, cluster_forming_p, df2, all_clusters, num_bins, this_sign);
                    this_prop = regionprops(CC,"Area");
                    size_distribution = [size_distribution; [this_prop.Area]'];
                else
                    [~, ~, ~, ~, CC, ~] = D2_cluster_detection_4D(this_dir, cluster_forming_p, df2, all_clusters, num_bins, this_sign);
                    this_prop = regionprops(CC,"Area");
                    this_prop_areas = [this_prop.Area];
                    this_prop_sorted = sort(this_prop_areas, "descend");
                    size_distribution = [size_distribution; this_prop_sorted(1:per_cluster)];
                end
            end
            total_number = size(size_distribution,1);
            [size_distribution_sorted, ~] = sort(size_distribution,1); % so that B = A(I)
            cutoff = size_distribution_sorted(round(total_number*(1-perm_p_val)),:);
            disp(['at p-value of ' num2str(perm_p_val) ' the smallest significant cluster has a size of ' num2str(cutoff(1))])
            save(['distribution_comparison_ttests_negative_s4_' num2str(total_number) 'perms'],'size_distribution_sorted');

            %% extract clusters
        case 3

            if ~isempty(load_dist)
                cd(perm_dir);
                load(load_dist);
                if redistribute_dist == 1
                    size_distribution_sorted = sort(reshape(size_distribution_sorted, 1, []))';
                elseif redistribute_dist == 2
                    size_distribution_sorted = size_distribution_sorted(:,1);
                end
                perm_num = size(size_distribution_sorted,1);
                cutoff = size_distribution_sorted(round(perm_num*(1-perm_p_val)),:);
                disp(['at p-value of ' num2str(perm_p_val) ' the smallest significant cluster has a size of ' num2str(cutoff)])
            end

            disp('applying cutoff to results of analysis...')
            if per_cluster == 0
                [realmaxArea, keepers, L, file_4D, CC, info] = D2_cluster_detection_4D(aoi_dir, cluster_forming_p, df2, cutoff, num_bins, this_sign);
            elseif save_3D_masks == 1 && redistribute_dist == 0
                [realmaxArea, keepers, L, file_4D, CC, info] = D2_cluster_detection_4D(aoi_dir, cluster_forming_p, df2, cutoff(3), num_bins, this_sign);
            else
                [realmaxArea, keepers, L, file_4D, CC, info] = D2_cluster_detection_4D(aoi_dir, cluster_forming_p, df2, 1, num_bins, this_sign);
            end

            if plot_result == 1 && (per_cluster == 0 || redistribute_dist > 0)

                plot([size_distribution_sorted]);
                hold on
                xline(round(perm_num*(1-perm_p_val)), 'red')
                if exp_plot == 1
                    set(gca, 'YScale', 'log')
                end
                h1 = get(gca, 'Children');
                legend([h1(1), h1(2)], {['alpha-level cut-off ' num2str(perm_p_val)], 'size of biggest cluster per permutation'})

                figure()
                plot([size_distribution_sorted; realmaxArea], 'o');
                hold on
                set(gca, 'YScale', 'log')
                yline(cutoff, 'red')
                h2 = get(gca, 'Children');
                legend([h2(1), h2(2)], {['minimal cluster size according to alpha-level cut-off: ' num2str(cutoff)], 'size of maximum cluster per analysis'})

            elseif plot_result == 1

                plot([size_distribution_sorted]);
                hold on
                xline(round(perm_num*(1-perm_p_val)), 'red')
                if exp_plot == 1
                    set(gca, 'YScale', 'log')
                end
                h1 = get(gca, 'Children');
                legend([h1(1), h1(2)], {['alpha-level cut-off ' num2str(perm_p_val)], ['size of ' num2str(per_cluster) ' biggest clusters per permutation']})

                prop = regionprops(CC,"Area");
                [cluster_sizes, sortingIdx] = sort([prop.Area], 'descend'); % so that B = A(I)

                for c = 1:per_cluster

                    thismaxArea = cluster_sizes(c);
                    figure()
                    plot([size_distribution_sorted(:,c); thismaxArea], 'o');
                    hold on
                    set(gca, 'YScale', 'log')
                    yline(cutoff(c), 'red')
                    h2 = get(gca, 'Children');
                    legend([h2(1), h2(2)], {['minimal cluster size according to alpha-level cut-off: ' num2str(cutoff(c))], ['size of ' num2str(c) '. biggest cluster per analysis: ' num2str(thismaxArea)]})

                end

            end

            if save_3D_masks

                cd(mask_save_dir);

                for b = 1:num_bins

                    this_3D_mask = cast(keepers(:,:,:,b),'single');
                    niftiwrite(this_3D_mask, ['mask_3D_3clusters_bin-' num2str(b)], info);

                end
            end

            %% print info about AOI
        case 4

            prop = regionprops(CC,"Area");
            [cluster_sizes, sortingIdx] = sort([prop.Area], 'descend'); % so that B = A(I)
            if per_cluster == 0
                kickIdx = sortingIdx(cluster_sizes<cutoff);
                idxs = setdiff(sortingIdx,kickIdx, 'stable');
            else
                kickIdx = sortingIdx([cluster_sizes(1:per_cluster)<cutoff cluster_sizes((per_cluster+1):end)<cutoff(end)]);
                idxs = setdiff(sortingIdx,kickIdx, 'stable');
            end

            if within_cluster_stats == 0

                to_visualize = zeros(size(file_4D));
                cluster_size = [];
                cluster_p_val = [];
                X_Y_Z = [];
                first_time = [];
                peak_time = [];
                last_time = [];
                peak_bin = [];
                peak_t = [];
                peak_p = [];
                counter = 0;
                regions = {};

                for i = idxs

                    counter = counter + 1;
                    otherwise_empty_brain = zeros(size(file_4D));

                    %           - size
                    %           - cluster p-value (based on comparison distribution)
                    %           - peak voxel (x y z): voxel and
                    %           - peak time point (bin / s)
                    %           - peak F-value

                    cluster_size = [cluster_size; prop(i).Area];
                    if per_cluster == 0 || redistribute_dist > 0
                        cluster_p_val = [cluster_p_val; 1-(find(size_distribution_sorted<prop(i).Area,1,'last')/(perm_num))];
                    else
                        cluster_p_val = [cluster_p_val; 1-(find(size_distribution_sorted(:,counter)<prop(i).Area,1,'last')/(perm_num))];
                    end

                    otherwise_empty_brain(L==i) = file_4D(L==i);
                    to_visualize(L==i) = file_4D(L==i);
                    [peakVal, peakVox] = max(otherwise_empty_brain, [], 'all');
                    [X, Y, Z, T] = ind2sub(size(file_4D), peakVox);
                    coordsMNI = ([X-1, Y-1, Z-1, 1] * info.Transform.T);
                    this_region = get_region_from_anatomy(coordsMNI(1:3)', 1);
                    [~,~,~,t1] = ind2sub(size(file_4D), find(otherwise_empty_brain,1,'first'))
                    [~,~,~,t2] = ind2sub(size(file_4D), find(otherwise_empty_brain,1,'last'))

                    X_Y_Z = [X_Y_Z; coordsMNI(1:3)]; % spm starts indexing at 0, apparently
                    first_time = [first_time; (t1-1)*0.18];
                    peak_time = [peak_time; (T-1)*0.18];
                    last_time = [last_time; (t2-1)*0.18];
                    peak_bin = [peak_bin; T];
                    peak_t = [peak_t; peakVal];
                    peak_p = [peak_p; 1-tcdf(peakVal,df2)];
                    regions = {regions{:}, this_region{:}};

                end

                peak_regions = regions';
                result_table = table(cluster_size, cluster_p_val, X_Y_Z, first_time, peak_time, last_time, peak_bin, peak_t, peak_p, peak_regions)

            else

                to_visualize = zeros(size(file_4D));
                cluster_size_4D = [];
                cluster_p_val = [];
                X_Y_Z_4D = [];
                peak_time = [];
                peak_bin = [];
                peak_t_4D = [];
                peak_p_4D = [];
                counter = 0;

                for i = idxs

                    counter = counter + 1;

                    disp(['4D Cluster ' num2str(counter) ':'])

                    otherwise_empty_brain = zeros(size(file_4D));

                    %           - size
                    %           - cluster p-value (based on comparison distribution)
                    %           - peak voxel (x y z): voxel and
                    %           - peak time point (bin / s)
                    %           - peak F-value

                    cluster_size_4D = [cluster_size_4D; prop(i).Area];
                    if per_cluster == 0
                        cluster_p_val = [cluster_p_val; 1-(find(size_distribution_sorted(:,counter)<prop(i).Area,1,'last')/(perm_num))];
                    else
                        cluster_p_val = [cluster_p_val; 1-(find(size_distribution_sorted(:,counter)<prop(i).Area,1,'last')/(perm_num))];
                    end

                    mask_4D = zeros(size(file_4D));
                    mask_4D(L == i) = 1;

                    otherwise_empty_brain(L==i) = file_4D(L==i);
                    to_visualize(L==i) = file_4D(L==i);
                    [peakVal, peakVox] = max(otherwise_empty_brain, [], 'all');
                    [X, Y, Z, T] = ind2sub(size(file_4D), peakVox);
                    coordsMNI = ([X-1, Y-1, Z-1, 1] * info.Transform.T);
                    this_region = get_region_from_anatomy(coordsMNI(1:3)', 1);

                    X_Y_Z_4D = [X_Y_Z_4D; coordsMNI(1:3)]; % spm starts indexing at 0, apparently
                    peak_time = [peak_time; (T-1)*0.18];
                    peak_bin = [peak_bin; T];
                    peak_t_4D = [peak_t_4D; peakVal];
                    peak_p_4D = [peak_p_4D; 1-tcdf(peakVal,df2)];


                    for b = 1:num_bins

                        disp(['3D Clusters in bin ' num2str(b) ':'])

                        this_3D = otherwise_empty_brain(:,:,:,b);
                        this_3D_mask = mask_4D(:,:,:,b);

                        [global_peakVal, global_peakVox] = max(this_3D, [], 'all');
                        [X, Y, Z] = ind2sub(size(this_3D), global_peakVox);
                        coordsMNI = ([X-1, Y-1, Z-1, 1] * info.Transform.T);
                        this_region = get_region_from_anatomy(coordsMNI(1:3)', 1);
                        disp(['Global Peak: ' num2str(coordsMNI(1:3)) ' (' this_region{:} ') with F-value ' num2str(global_peakVal)])


                        CC_3D = bwconncomp(this_3D_mask);
                        L_3D = labelmatrix(CC_3D);
                        prop_3D = regionprops(CC_3D,"Area");
                        [cluster_sizes_3D, sortingIdx_3D] = sort([prop_3D.Area], 'descend'); % so that B = A(I)

                        if ~isempty(min_k_3D)
                            idx_3D = sortingIdx_3D(cluster_sizes_3D>min_k_3D);
                        else
                            idx_3D = sortingIdx_3D;
                        end

                        cluster_size_3D = []; % cluster_size_3D(sortingIdx_3D(idx_3D));
                        X_Y_Z_3D = [];
                        peak_t_3D = [];
                        peak_p_3D = [];
                        region_3D = {};

                        for c3d = idx_3D

                            cluster_size_3D = [cluster_size_3D; prop_3D(c3d).Area];
                            otherwise_empty_3D = zeros(size(this_3D));

                            otherwise_empty_3D(L_3D == c3d) = this_3D(L_3D == c3d);
                            this_cluster = this_3D(L_3D == c3d);

                            [cluster_peakVal, cluster_peakVox] = max(otherwise_empty_3D, [], 'all');
                            [X, Y, Z] = ind2sub(size(otherwise_empty_3D), cluster_peakVox);
                            coordsMNI = ([X-1, Y-1, Z-1, 1] * info.Transform.T);
                            this_region = get_region_from_anatomy(coordsMNI(1:3)', 1);


                            X_Y_Z_3D = [X_Y_Z_3D; coordsMNI(1:3)]; % spm starts indexing at 0, apparently
                            peak_t_3D = [peak_t_3D; cluster_peakVal];
                            peak_p_3D = [peak_p_3D; 1-tcdf(cluster_peakVal,df2)];
                            region_3D = {region_3D{:}, this_region{:}};

                        end

                        peak_region3D = region_3D';
                        result_table = table(cluster_size_3D, X_Y_Z_3D, peak_t_3D, peak_p_3D, peak_region3D)
                        answer = input('More?');

                        if answer == 0 % I've seen enough

                            break;

                        end

                    end

                end

                result_table = table(cluster_size_4D, cluster_p_val, X_Y_Z_4D, peak_time, peak_bin, peak_t_4D, peak_p_4D)

            end
    end

end
