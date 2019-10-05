function [ Pr] = Channel_with_Memory (MODE  , epsilon_1 , epsilon_2 , numLevel , noise_index , delta)
%% Finding the channel transition probaility matrix Pr(j + 1 | i + 1 )

Pr_z_1i_given_z_1i_1 = [(1 - epsilon_1(noise_index) + delta) / (1 + delta)  , epsilon_1(noise_index) / (1 + delta) ;
    (1 - epsilon_1(noise_index)) / (1 + delta) , (epsilon_1(noise_index) + delta ) / (1 + delta)] ;

Pr_z_2i_given_z_2i_1 = [(1 - epsilon_2(noise_index) + delta) / (1 + delta)  , epsilon_2(noise_index) / (1 + delta) ;
    (1 - epsilon_2(noise_index)) / (1 + delta) , (epsilon_2(noise_index) + delta ) / (1 + delta)] ;

Pr = zeros(numLevel ^ 2 , numLevel ^ 2) ;
%% Alternative
n = log2(numLevel)  ;
for x_1 = 1 : numLevel
    for x_2 = 1 : numLevel
        for y_1 = 1 : numLevel
            for y_2 = 1 : numLevel
                binary_x_1 = de2bi(x_1 - 1 , n , 'left-msb') ;
                binary_x_2 = de2bi(x_2 - 1 , n , 'left-msb') ;
                binary_y_1 = de2bi(y_1 - 1 , n , 'left-msb') ;
                binary_y_2 = de2bi(y_2 - 1 , n , 'left-msb') ;
                if (MODE == 1)
                    binary_z_1 = xor(binary_x_2 , binary_y_1) ;
                    binary_z_2 = xor(binary_x_1 , binary_y_2) ;
                elseif (MODE == 2)
                    hold_var = xor(binary_x_1 , binary_x_2) ;
                    binary_z_1 = xor(hold_var , binary_y_1) ;
                    binary_z_2 = xor(hold_var , binary_y_2) ;
                elseif (MODE == 3)
                    hold_var = binary_x_1 .* binary_x_2 ;
                    binary_z_1 = xor(hold_var , binary_y_1) ;
                    binary_z_2 = xor(hold_var , binary_y_2) ;
                end
                if (binary_z_1(1) == 0 && binary_z_2(1) == 0)
                    myPr = (1 - epsilon_1(noise_index)) * (1 - epsilon_2(noise_index)) ;
                    
                elseif (binary_z_1(1) == 0 && binary_z_2(1) == 1)
                    myPr = (1 - epsilon_1(noise_index)) * epsilon_2(noise_index) ;
                    
                elseif (binary_z_1(1) == 1 && binary_z_2(1) == 0)
                    myPr = epsilon_1(noise_index) * (1 - epsilon_2(noise_index)) ;
                    
                else
                    myPr = epsilon_1(noise_index) * epsilon_2(noise_index) ;
                end
                
                for i = 2 : n
                    myPr = myPr * Pr_z_1i_given_z_1i_1(binary_z_1(i - 1) + 1 , binary_z_1(i) + 1) ...
                        * Pr_z_2i_given_z_2i_1(binary_z_2(i - 1) + 1 , binary_z_2(i) + 1) ;
                end
                x = (x_2 - 1) * numLevel + x_1 ;
                y = (y_1 - 1) * numLevel + y_2 ;
                
                Pr(x , y) = myPr ;
            end
        end
    end
end
end