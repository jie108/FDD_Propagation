function [grid_prob] = Grid_prob(grid_u, dir_v_norm, k, thre)
%% parameters: FOD: single peak FOD: n by 1 vector 
 %% grid_u:  the corresponding grid points, 3 by n matrix; dir_v_norm: normalized direction vectors (26)
 %% k: sharing between the k -closest directions 
 %% thre: cutoff value to thresholding small prob. ; should be no larger than 1/26
 %% return: FDD: 26 by 1 vector 
%%return: FDD vector: normalized according to group size and sum to 1
grid_prob=zeros(size(grid_u,2),26);


 for i = 1: size(grid_u,2)
     grid_c=grid_u(:,i);
     inner_c=dir_v_norm*grid_c;
     dist_c=acos(inner_c);
     
     index_c= find(dist_c==min(dist_c));
     
     if min(dist_c)==0
         grid_prob(i,index_c)=1;
     else
         [temp I]=sort(dist_c); 
         grid_prob(i,I(1:k))=1./dist_c(I(1:k))./sum(1./dist_c(I(1:k)));
         grid_prob(i,grid_prob(i,:)<thre)=0;
     end
 end
%%% normalize by group size and sum to 1