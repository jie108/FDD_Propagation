function [group_id group_size dist_min] =dir_group(pos_plot, dir_v_norm)
%%parameter: pos_plot: 3d coordinates of the points on a grid 
%% dir_v_norm: unit-norm direction vectors
%% return: group_id: which direction vector each grid point is closest to 
%% group_size: size of each group
%% dist_min: minimum angular distance with the FDD directions


group_id=zeros(size(pos_plot,2),1);
dist_min=zeros(size(pos_plot,2),1);

for i=1:size(pos_plot,2)
    grid_c=pos_plot(:,i);
    ang_c=acos(dir_v_norm*grid_c);
    group_id(i,1)=find(ang_c==min(ang_c));
    dist_min(i,1)=min(ang_c);
end

group_size=zeros(size(dir_v_norm,1),1); %%size of each region

for i=1:size(group_id,1)
    for j = 1:size(dir_v_norm,1)
        if group_id(i,1)==j;
            group_size(j,1)=group_size(j,1)+1;
        end
    end
end