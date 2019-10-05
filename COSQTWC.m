function [SDR , SDR_1 , SDR_2 , joint_T , codebook_user_1 , codebook_user_2 , Distortion] = COSQTWC ( p , Pr , numLevel , joint_T , width , f , f_u_2_given_u_1 , f_u_1_given_u_2)
delta_A = width (2) * width (1) ;
codebook_user_1 = zeros (numLevel , 1 , length(joint_T)) ;
codebook_user_2 = zeros (numLevel , 1 , length(joint_T)) ;
Threshold = 0.001 ;
D = [1 2] ;
FileID = fopen ('Results.txt', 'a') ;

while abs(D(2) - D(1)) / D(2) > (Threshold)
    D (1) = D(2) ;
    %% Optimal centroids for user 1
    parfor u_2_index = 1 : length(joint_T)
        u_2 = joint_T (u_2_index , 3) ;
        for j_2 = 1 : numLevel
            numerator = 0 ;
            denominator = 0 ;
            for j_1 = 1 : numLevel
                for i_1 = 1 : numLevel
                    u_1_index = find (joint_T(: , 2) == i_1) ;
                    u_1 = joint_T(u_1_index , 1) ;
                    
                    i = (joint_T (u_2_index , 4) - 1) * numLevel + i_1 ;
                    j = (j_1 - 1) * numLevel + j_2 ;
                    
                    numerator = numerator + Pr ( i , j) * sum (u_1 .* f_u_1_given_u_2(u_1_index , u_2_index)) ;
                    denominator = denominator + Pr (i , j) * sum (f_u_1_given_u_2(u_1_index , u_2_index)) ;
                end
            end
            if (numerator == 0 && denominator == 0)
                
                codebook_user_1 (j_2 , : , u_2_index) = p * u_2 ;
            else
                codebook_user_1 (j_2 , : , u_2_index) = numerator / denominator ;
            end
        end
    end
    
    %% Optimal centroids for user 2
    parfor u_1_index = 1 : length(joint_T)
        u_1 = joint_T (u_1_index , 1) ;
        for j_1 = 1 : numLevel
            numerator = 0 ;
            denominator = 0 ;
            for j_2 = 1 : numLevel
                for i_2 = 1 : numLevel
                    u_2_index = find (joint_T (: , 4) == i_2) ;
                    u_2 = joint_T (u_2_index , 3) ;
                    
                    i = (i_2 - 1) * numLevel + joint_T (u_1_index , 2) ;
                    j = (j_1 - 1) * numLevel + j_2 ;
                    numerator = numerator + Pr (i , j) * sum (u_2 .* f_u_2_given_u_1(u_1_index , u_2_index)') ;
                    denominator = denominator + Pr (i , j) * sum (f_u_2_given_u_1(u_1_index , u_2_index)') ;
                end
            end
            if (numerator == 0 && denominator == 0)
                
                codebook_user_2 (j_1 , : , u_1_index) = p * u_1 ;
            else
                codebook_user_2 (j_1 , : , u_1_index) = numerator / denominator ;
            end
        end
    end
    
    %% Optimal Partitions for user 1
    joint_T_u_1 = zeros(length(joint_T) , 1) ;
    parfor u_1_index = 1 : length(joint_T)
        temp = zeros (1 , numLevel) ;
        summation = 0 ;
        u_1 = joint_T (u_1_index , 1) ;
        for i_1 = 1 : numLevel
            for j_1 = 1 : numLevel
                for j_2 = 1 : numLevel
                    for i_2 = 1 : numLevel
                        
                        i = (i_2 - 1) * numLevel + i_1 ;
                        j = (j_1 - 1) * numLevel + j_2 ;
                        u_2_index = find (joint_T (: , 4) == i_2) ;
                        u_2 = joint_T (u_2_index , 3) ;
                        
                        
                        hold_var = u_1 - codebook_user_1 (j_2 , : , u_2_index) ;
                        hold_var = hold_var (:) ;
                        
                        summation = summation + Pr (i , j) * width(2) ...
                            * sum (f_u_2_given_u_1(u_1_index , u_2_index)' .* ((hold_var) .^ 2 ...
                            + (u_2 - codebook_user_2 (j_1 , : , u_1_index)) .^ 2)) ;
                    end
                end
            end
            temp (i_1) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min (temp) ;
        joint_T_u_1 (u_1_index , 1) = partition_index ;
    end
    joint_T(: , 2) = joint_T_u_1 (: , 1) ;
    
    %% Optimal paartitions for user 2
    joint_T_u_2 = zeros(length(joint_T) , 1) ;
    parfor u_2_index = 1 : length (joint_T)
        temp = zeros (1 , numLevel) ;
        summation = 0 ;
        u_2 = joint_T (u_2_index , 3) ;
        for i_2 = 1 : numLevel
            for j_1 = 1 : numLevel
                for j_2 = 1 : numLevel
                    for i_1 = 1 : numLevel
                        i = (i_2 - 1) * numLevel + i_1 ;
                        j = (j_1 - 1) * numLevel + j_2 ;
                        
                        u_1_index = find (joint_T (: , 2) == i_1) ;
                        u_1 = joint_T (u_1_index , 1) ;
                        
                        
                        hold_var = u_2 - codebook_user_2 (j_1 , : , u_1_index) ;
                        hold_var = hold_var (:) ;
                        
                        summation = summation + Pr (i , j) * width (1) * sum (f_u_1_given_u_2(u_1_index , u_2_index) .* ((u_1 - codebook_user_1 (j_2 , : , u_2_index)) .^ 2 ...
                            + (hold_var) .^ 2)) ;
                    end
                end
            end
            temp (i_2) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min (temp) ;
        joint_T_u_2 (u_2_index , 1) = partition_index ;
    end
    
    joint_T(: , 4) = joint_T_u_2 (: , 1) ;
    [D(2) , D_1 , D_2] = distortion (Pr , joint_T , codebook_user_1 , codebook_user_2 , numLevel, delta_A , f) ;
    fprintf (FileID , '\n Overall: \n') ;
    fprintf (FileID , ' %7.4f ' , D(2)) ;
    fprintf (FileID , ' %7.4f ' , D_1 ) ;
    fprintf (FileID , ' %7.4f ' , D_2 ) ;
end
Distortion = D(2) ;
SDR = 10 * log10(2 / D(2)) ;
SDR_1 = 10 * log10(1 / D_1) ;
SDR_2 = 10 * log10(1 / D_2) ;
fclose (FileID) ;
end