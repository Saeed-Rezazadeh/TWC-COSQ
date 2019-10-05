% This is the MATLAB script of the proposed Algorithm 3 in the Thesis. 
clc;
close all
clear all
%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon. 
% Also, since the proposed ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for given a epsilon and delta. 
FileID = fopen ('Results.txt', 'a') ;

%% Number of Quantization Levels
QuantizationLevels = [4 16] ;


rho = [0 , 0.5 , 0.9 0.99] ;

%% Noise parameters
% Channel MODE
MODE = 3 ; % 1 = Orthogonal  2 = Additive 3 = Multiplicative

%% Channel's cross-over probability epsilon for eeach direction of transmission. 
epsilon_1 = [0  0.005 0.01 0.05] ; 
epsilon_2 = [0  0.01 0.05 0.1] ; 

%% The variable delta determines the amount of noise correlation. Choose between 0, 5, and 10
delta = 10 ;

% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method. 
SIZE = length(epsilon_1) ;
i_index = [1 : SIZE , SIZE : -1 : 1  , 1 : SIZE , SIZE : -1 : 1] ;
%% SDR values
SDR = zeros (length(QuantizationLevels) , length(rho) , length(i_index)) ;
final_SDR = zeros (length(QuantizationLevels) , length(rho) , SIZE) ;

SDR_1 = zeros (length(QuantizationLevels) , length(rho) , length(i_index)) ;
final_SDR_1 = zeros (length(QuantizationLevels) , length(rho) , SIZE) ;

SDR_2 = zeros (length(QuantizationLevels) , length(rho) , length(i_index)) ;
final_SDR_2 = zeros (length(QuantizationLevels) , length(rho) , SIZE) ;


for l = 1 : length (QuantizationLevels)
    numLevel = QuantizationLevels (l);
    for j = 1  : length(rho)
        p = rho (j) ;
        u = [0 0] ;
        for k = 1 : length(i_index)
            i = i_index(k) ;
            Pr = Channel_with_Memory(MODE , epsilon_1 , epsilon_2 , numLevel , i , delta) ;
            
            %% Uncomment this section if the variable MODE is set to 3. 
            % Initial partitions obtained from the half-duplex design. 
            if (k == 1)
              % the following initialization is particularly for BM-TWC
              % whose details are provided in the Thesis Section. 3.2.2. 
              if (numLevel == 16)
                  clear joint_T
                  Data = ['BM-TWC\delta_' num2str(delta) '\Data_' num2str(k) '_p_' num2str(j) '_numLevel_' num2str(4)];
              else
                  clear joint_T
                  Data = ['BM-TWC\delta_' num2str(delta) '\Data_' num2str(k) '_p_' num2str(j) '_numLevel_' num2str(2)];
              end
              load (Data) ;
            end
            
            %% Uncomment this section if the variable MODE is set to 2. 
            % Initial partitions obtained from the half-duplex design. 
            % Data = ['BA-TWC\delta_' num2str(delta) '\Data_' num2str(k) '_p_' num2str(j) '_numLevel_' num2str(numLevel)];
            % clear joint_T ; 
            % load (Data) ; 
            
            
            width (2) = joint_T (2 , 3) - joint_T (1 , 3) ;
            width (1) = joint_T (2 , 1) - joint_T (1 , 1) ;
            
            [f , f_u_2_given_u_1 , f_u_1_given_u_2] = generate_PDF(joint_T , p , width) ;
            
            
            [SDR(l , j , k)  , SDR_1(l , j , k)  , SDR_2(l , j , k)  , joint_T , codebook_user_1 , codebook_user_2] ...
                = COSQTWC ( p , Pr , numLevel , joint_T , width , f , f_u_2_given_u_1 , f_u_1_given_u_2) ;
            
            fprintf (FileID , '\n i = %4.2f \n' , i ) ;
            fprintf (FileID , '\n SDR values: \n') ;
            fprintf (FileID , ' SDR = %7.4f SDR1 = %7.4f SDR2 = %7.4f' , SDR(l , j , k) , SDR_1(l , j , k) , SDR_2(l , j , k)) ;
            Data = ['T_' num2str(delta) '\Main_Data_' num2str(k) '_p_' num2str(j) '_numLevel_' num2str(numLevel) '_delta_' num2str(delta)] ;
            save(Data , 'joint_T' , 'codebook_user_1' , 'codebook_user_2') ;
        end
        
        fprintf (FileID , '\nOver All SDR\n') ;
        for i = 1 : SIZE
            index = find (i_index == i) ;
            hold_var = SDR(l , j , index) ;
            hold_var = hold_var(:) ;
            final_SDR(l , j , i) = max(hold_var) ;
            fprintf (FileID , '\ni = %d' , i) ;
            fprintf (FileID , '\nSDR = %7.4f' , final_SDR(l , j , i)) ;
        end
        fprintf (FileID , '\nSDR for user 1\n') ;
        for i = 1 : SIZE
            index = find (i_index == i) ;
            hold_var = SDR_1(l , j , index) ;
            hold_var = hold_var(:) ;
            final_SDR_1(l , j , i) = max(hold_var) ;
            fprintf (FileID , '\nSDR_1 = %7.4f\n' , final_SDR_1(l , j , i)) ;
        end
        fprintf (FileID , '\nSDR for user 2\n') ;
        for i = 1 : SIZE
            index = find (i_index == i) ;
            hold_var = SDR_2(l , j , index) ;
            hold_var = hold_var(:) ;
            final_SDR_2(l , j , i) = max(hold_var) ;
            fprintf (FileID , '\nSDR_2 = %7.4f\n' , final_SDR_2(l , j , i)) ;
        end
        
        Data = ['Main_Data' '_p_' num2str(j) '_numLevel_' num2str(numLevel)] ;
        save(Data , 'SDR' , 'SDR_1' , 'SDR_2') ;
        
        fprintf (FileID , '\n p = %4.2f \n' , p ) ;
    end
    fprintf (FileID , '\n numLevel = %d \n' , numLevel ) ;
end

