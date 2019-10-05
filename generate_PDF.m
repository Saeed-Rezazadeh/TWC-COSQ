function  [f , f_u_2_given_u_1 , f_u_1_given_u_2] = generate_PDF (joint_T , p , width)
f_u_2_given_u_1 = zeros(length(joint_T) , length(joint_T)) ; 
f = zeros(length(joint_T) , length(joint_T)) ; 
u_2 = joint_T(: , 3) ;
for u_1_index = 1 : length(joint_T)
    u_1 = joint_T(u_1_index , 1) ;
    
    f_u_2_given_u_1(u_1_index , :) =  1 ./ (sqrt (2 .* pi .* (1 - p .^ 2))) .* exp (- 1 ./ (2 .* (1 - p .^ 2)) .* (u_2 - p .* u_1) .^ 2) ;
    f(u_1_index , :) = 1 ./ (2 .* pi .* sqrt (1 - p .^ 2)) .* exp (-1 ./ (2 .* (1 - p .^ 2)) .* (u_1 .^ 2 + u_2 .^ 2 - 2 .* p .* u_1 .* u_2)) ;
    
end

f_u_1_given_u_2 = zeros(length(joint_T) , length(joint_T)) ; 
u_1 = joint_T(: , 1) ;
for u_2_index = 1 : length(joint_T) 
    u_2 = joint_T(u_2_index , 3) ; 
    f_u_1_given_u_2(: , u_2_index) = 1 ./ (sqrt (2 .* pi .* (1 - p .^ 2))) .* exp (- 1 ./ (2 .* (1 - p .^ 2)) .* (u_1 - p .* u_2) .^ 2) ;
end 

f_u_1_given_u_2 = f_u_1_given_u_2 ./ (sum(f_u_1_given_u_2 , 1) * width (1)) ; 
f_u_2_given_u_1 = f_u_2_given_u_1 ./ (sum(f_u_2_given_u_1 , 2) * width (2)) ; 
f = f./ (sum(sum(f) .* width(1)) * width(2)) ; 
end