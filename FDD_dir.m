function [dir_v dir_v_norm]= FDD_dir()
%% return: 26 FDD directions
dir_v=zeros(26,3);  %%26 FDD directions
dir_v_norm=dir_v;   %% normalized FDD directions
count=1;
for i=-1:1
    for j=-1:1
        for k=-1:1
            cur_v=[i j k];
            cur_n=sqrt(sum(cur_v.^2));
            if abs(cur_n)>1e-6
                dir_v(count,:)= cur_v;
                dir_v_norm(count,:)=cur_v./cur_n;
                count=count+1;
            end
        end
    end
end
