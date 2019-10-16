function [FDD] = FOD_to_FDD(FOD, group_id, group_size)
%% parameters: FOD: FOD evaluated on a grid, nonnegative 
%% group_id: which group each grid points belongs to
%% group_size: size of each group
%%return: FDD vector: normalized according to group size and sum to 1

FDD=zeros(size(group_size,1),1);


 for i = 1: size(group_id,1)
     for j = 1: size(group_size,1)
         if group_id(i,1)==j
             FDD(j,1)=FDD(j,1)+FOD(i,1);
         end
     end
 end
%%% normalize by group size and sum to 1
FDD=FDD./group_size;
FDD=FDD./sum(FDD);