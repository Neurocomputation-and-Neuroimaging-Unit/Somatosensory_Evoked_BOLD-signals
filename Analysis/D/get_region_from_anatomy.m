% the automatic region retreiver using anatomy toolbox
% input: x y z in MNI, top N results
% output: top N regions + according probability

function [top_probN] = get_region_from_anatomy(xyzmm, n)

addpath('...\Anatomy_v3')
load('...\Anatomy_v3\JuBrain_Data_public_v30.mat');
xyz = JuBrain.Vo.mat \ [xyzmm; 1];
xyz = round(xyz(1:3));

[loc, ~] = se_CytoForVoxel(xyz,JuBrain,5);

try 
    top_probN = {loc{2,1}{1:n,:}};
catch
    top_probN = {'unassigned'};
end
end

function [loc, ref] = se_CytoForVoxel(xyzv,Atlas,cut)
global st

index = Atlas.Index(min(max(xyzv(1),1),size(Atlas.Index,1)),min(max(xyzv(2),1),size(Atlas.Index,2)),min(max(xyzv(3),1),size(Atlas.Index,3)));

ref = {'',''};

if index == 0
    loc = {'',''};
    return
    
else
    
    try
        idx = Atlas.mpm(index);
        loc{1,1} = Atlas.Namen{idx};
        ref1 = Atlas.Paper{idx};
        
    catch
        loc{1,1} = '';
        ref1 = [];
    end
    
    ref2 = [];
    
    [a, ~, B] = find(Atlas.PMap(:,index)); [B I] = sort(B,'descend'); B = B(B>.01);
    if numel(B)>0
        for i=1:min(cut,numel(B))
            tmp{i,1} = [sprintf('%4.1f',B(i)*100) '%   ' Atlas.Namen{a(I(i))}];
            if isfield(Atlas,'Paper')
                ref2 = [ref2 Atlas.Namen{a(I(i))} '\n'];
                ref2 = [ref2 Atlas.Paper{a(I(i))} '\n'];
            end
        end
        loc{2,1} = tmp;
    else
        loc{2,1} = '';
    end
    ref = {ref1, ref2};
end
end
