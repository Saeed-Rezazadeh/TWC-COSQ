function [ width , joint_T ] = generate_source ( u , p , alpha , resolution , numLevel)
Var = [1 p ;
    p 1] ;
R = mvnrnd(u,Var,alpha);

T_1 = R(: , 1) ; 
T_2 = R(: , 2) ; 

MIN_1 = min(T_1) ;
MAX_1 = max(T_1) ;


MIN_2 = min(T_2) ;
MAX_2 = max(T_2) ;

width(1) = (MAX_1 - MIN_1) / resolution ;
width(2) = (MAX_2 - MIN_2) / resolution ;

u_1_edge = (MIN_1 : width(1) : MAX_1)' ;
u_2_edge = (MIN_2 : width(2) : MAX_2)' ;

initial_codebook_user_1 = datasample(T_1 , numLevel ) ;
[partition_user_1 , ~ , ~] = lloyds(T_1 ,initial_codebook_user_1);


u_1_edge (: , 2) = quantiz (u_1_edge , partition_user_1) + 1 ;

initial_codebook_user_2 = datasample(T_2 , numLevel ) ;
[partition_user_2 , ~ , ~] = lloyds(T_2 ,initial_codebook_user_2);

u_2_edge (: , 2) = quantiz (u_2_edge , partition_user_2) + 1; 
joint_T = cat (2 , u_1_edge , u_2_edge) ; 
[~ , ~ , ~ , joint_T , ~ , ~] = COSQTWC ( p , eye (numLevel ^ 2 , numLevel ^ 2) , numLevel , joint_T , width ) ; 
end
