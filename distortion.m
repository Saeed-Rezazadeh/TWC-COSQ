function [D , D_1 , D_2] = distortion (Pr , joint_T , codebook_user_1 , codebook_user_2 , numLevel, delta_A , f)
summation = 0 ;
summation_1 = 0 ;
summation_2 = 0 ;

parfor i_1 = 1 : numLevel
    for i_2 = 1 : numLevel
        for j_1 = 1 : numLevel
            for j_2 = 1 : numLevel
                
                i = (i_2 - 1) * numLevel + i_1 ;
                j = (j_1 - 1) * numLevel + j_2 ;
                
                u_1_index = find (joint_T (: , 2) == i_1) ;
                
                u_2_index = find (joint_T (: , 4) == i_2) ;
                u_2 = joint_T (u_2_index , 3) ;
                
                for u_1_i = 1 : length(u_1_index)
                    u_1 = joint_T (u_1_index (u_1_i) , 1) ;
                    
                    hold_var = u_1 - codebook_user_1 (j_2 , : , u_2_index) ;
                    hold_var = hold_var (:) ;
                    
                    summation = summation + Pr (i , j) * delta_A ...
                        * sum (f(u_1_index(u_1_i) , u_2_index)' .* ((hold_var) .^ 2 ...
                        + (u_2 - codebook_user_2 (j_1 , : , u_1_index (u_1_i))) .^ 2)) ;
                    summation_1 = summation_1 + Pr (i , j) * delta_A * sum (f(u_1_index(u_1_i) , u_2_index)' .* (hold_var) .^ 2) ;
                    summation_2 = summation_2 + Pr (i , j) * delta_A * sum (f(u_1_index(u_1_i) , u_2_index)' .* (u_2 - codebook_user_2 (j_1 , : , u_1_index (u_1_i))) .^ 2) ;
                end
            end
        end
    end
end
D = summation ;
D_1 = summation_1 ;
D_2 = summation_2 ;
end