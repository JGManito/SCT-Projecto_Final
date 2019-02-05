%--------------------------------------------------------%
%           SCT - Projecto Final                         %
%--------------------------------------------------------%
%           Projecto 1                                   %
%                                                        %
%Pedro Afonso                                            %
%Joao Manito                                             %
%--------------------------------------------------------%
%Funcao main                                             %
%--------------------------------------------------------%

clear;
clc;

%Define problem variables
mask_angle = 10;
noise_sigma = sqrt(10);
N_PDOP = 5;
disable_canyon = 1;
disable_noise = 0;
disable_ionosphere = 1;
NSimulations = 100;

%Load ephemeris data
%Data used: Yuma week 1014 - 319488
yuma_dataset=read_yuma("almanac.yuma.week1014.319488.txt");

for iteration = 1:NSimulations
    %Initialize position
    pos_real_current_WGS84=[40,-9,2000];
    pos_real_current_ECEF=LLH2XYZ(pos_real_current_WGS84);
    
    %Initialize position and velocity in XY plane
    pos_real_current_xy=[0,0];
    velocity_real_current_xy=[0,0];
    
    %Define the initial position estimate as the center of Portugal
    position_estimate = LLH2XYZ([39.5, -8.0,0]);
    
    %Start the simulation
    for t_sim=0:1:250 %Total trajectory time is 250s, with 1Hz sample rate
        
        fprintf("Iteration #%d, time =%d\n",iteration,t_sim);
        
        %Change to N=4 for canyon scenario
        if t_sim >= 201 && disable_canyon == 0
            flag = 1;
        else
            flag = 0;
        end
        
        
        %------Aircraft Motion Simulation------%
        %Save previous position
        if t_sim == 0
            pos_real_prev_xy = pos_real_current_xy(1,:);
            pos_real_prev_WGS84 = pos_real_current_WGS84(1,:);
            pos_real_prev_ECEF = LLH2XYZ(pos_real_current_WGS84(1,:));
            pos_estimate_prev_ECEF = position_estimate;
        else
            pos_real_prev_xy = pos_real_current_xy(t_sim,:);
            pos_real_prev_WGS84 = pos_real_current_WGS84(t_sim,:);
            pos_real_prev_ECEF = LLH2XYZ(pos_real_current_WGS84(t_sim,:));
            pos_estimate_prev_ECEF = position_estimate_ECEF(t_sim,:);
        end
        
        %Update position in the XY plane
        [pos_real_current_xy(t_sim+1,:),velocity_real_current_xy(t_sim+1,:)]=update_position(t_sim);
        
        %Compute the new position in WGS84 and ECEF frames
        pos_real_current_WGS84(t_sim+1,:) = flat2llh(pos_real_current_xy(t_sim+1,:) , pos_real_prev_xy , pos_real_prev_WGS84);
        pos_real_current_ECEF(t_sim+1,:) = LLH2XYZ(pos_real_current_WGS84(t_sim+1,:));
        
        
        
        
        %------Satellite Motion Simulation------%
        %Get satellite positions in ECEF
        satellite = compute_orbital_parameters(yuma_dataset);
        [satellite_positions,satellite] = compute_satellite_position(satellite,t_sim);
        %disp(satellite_positions);
        
        %Compute the true range to the receiver
        true_range = compute_true_range(satellite_positions,pos_real_current_ECEF(t_sim+1,:));
        
        %Determine satellites in view using the initial estimate for position
        [satellite_positions_visible,satellite_visible,true_range_visible,elevation,N_view] = filter_satellites_visible_LS(satellite_positions,satellite,true_range,mask_angle,position_estimate,flag);
        NView(t_sim+1) = N_view;
        
        %Use PDOP
        [satellite_positions_minimized] = minimize_PDOP(satellite_positions_visible,position_estimate,N_PDOP);
        
        %Compute the pseudorange to the receiver
        pseudorange_visible = compute_pseudoranges(true_range_visible,elevation,noise_sigma,disable_ionosphere,disable_noise);
        
        
        
        %------Receiver Motion Simulation------%
        
        
        %Find the Least Squares Solution
        [position_estimate] = compute_LS_solution(satellite_positions_visible,position_estimate,pseudorange_visible);
        
        position_estimate_ECEF(t_sim+1,:) = position_estimate (1:3);
        position_estimate_WGS84(t_sim+1,:) = XYZ2LLH(position_estimate(1:3));
        
        
        N_XY(t_sim+1)=pos_real_current_xy(t_sim+1,2);
        E_XY(t_sim+1)=pos_real_current_xy(t_sim+1,1);
        
        N_LLH(t_sim+1)=pos_real_current_WGS84(t_sim+1,1);
        E_LLH(t_sim+1)=pos_real_current_WGS84(t_sim+1,2);
        
        N_estimate_LLH(t_sim+1)=position_estimate_WGS84(t_sim+1,1);
        E_estimate_LLH(t_sim+1)=position_estimate_WGS84(t_sim+1,2);
        
        
        VX(t_sim+1) = position_estimate_ECEF(t_sim+1,1) - pos_estimate_prev_ECEF(1);
        VY(t_sim+1) = position_estimate_ECEF(t_sim+1,2) - pos_estimate_prev_ECEF(2);
        VZ(t_sim+1) = position_estimate_ECEF(t_sim+1,3) - pos_estimate_prev_ECEF(3);
        
        VABS(t_sim+1) = sqrt(VX(t_sim+1)^2 + VY(t_sim+1)^2 + VZ(t_sim+1)^2);
        
        
    end
    
    %Error analysis
    %Absolute position error
    for i=1:251
        pos_real_current_ECEF(i,:) = LLH2XYZ(pos_real_current_WGS84(i,:));
        errorX(i) = pos_real_current_ECEF(i,1) - position_estimate_ECEF(i,1);
        errorY(i) = pos_real_current_ECEF(i,2) - position_estimate_ECEF(i,2);
        errorZ(i) = pos_real_current_ECEF(i,3) - position_estimate_ECEF(i,3);
        error_position(i) = sqrt(errorX(i)^2 + errorY(i)^2 + errorZ(i)^2);
    end
    
    %Absolute velocity error
    for i=1:251
        if i<101
            v_real = 50 + (i-1);
        else
            v_real = 150;
        end
        
        error_velocity(i) = VABS(i) - v_real;
    end
    
    %RMS position and velocity error
    error_position_RMS(iteration) = sqrt((sum(error_position.^2))/250);
    error_velocity_RMS(iteration) = sqrt((sum(error_velocity.^2))/250);
    %Compute the RMS error ignoring the initial estimation
    error_position_RMS_offset(iteration) = sqrt((sum(error_position(2:end).^2))/249);
    error_velocity_RMS_offset(iteration) = sqrt((sum(error_velocity(2:end).^2))/249);
    
    
    if NSimulations == 1
        
        plot(E_XY,N_XY);
        title('Valores reais plano NE');
        grid on
        daspect([1 1 1])
        
        figure
        plot(E_LLH,N_LLH);
        title('Valores reais LLH');
        grid on
        daspect([1 1 1])
        
        figure
        plot(E_estimate_LLH,N_estimate_LLH);
        title('Valores estimados LLH');
        grid on
        daspect([1 1 1])
        
        figure
        plot(0:1:250,VABS)
        title('Velocidade estimada');
        grid on
        
        figure
        semilogy(0:1:250,abs(error_position))
        title('Erro de posição');
        grid on
        
        figure
        semilogy(0:1:250,abs(error_velocity))
        title('Erro de velocidade');
        grid on
        
        figure
        plot(0:1:250,NView)
        title('Número de satélites em vista');
        grid on
        
        fprintf("The simulation had the following conditions:\n");
        if disable_canyon == 1
            fprintf("Canyon scenario disabled\n");
        else
            fprintf("Canyon scenario enabled\n");
        end
        if disable_ionosphere == 1
            fprintf("Ionospheric perturbations disabled\n");
        else
            fprintf("Ionospheric perturbations enabled\n");
        end
        
        fprintf("The noise sigma is: %f\nThe number of satellites to use for the PDOP minimization is %d\n",noise_sigma, N_PDOP);
        
        
        fprintf("The position RMS error for 1 iteration is %f\n",error_position_RMS);
        fprintf("The velocity RMS error for 1 iteration is %f\n",error_velocity_RMS);
        fprintf("The position RMS error for 1 iteration (ignoring the first position estimation) is %f\n",error_position_RMS_offset);
        fprintf("The velocity RMS error for 1 iteration (ignoring the first velocity estimation) is %f\n",error_velocity_RMS_offset);
    end
    
    
end

if NSimulations ~= 1
    
    %Compute the average RMS errors
    error_position_RMS_avg = sum(error_position_RMS)/NSimulations;
    error_velocity_RMS_avg = sum(error_velocity_RMS)/NSimulations;
    error_position_RMS_offset_avg = sum(error_position_RMS_offset)/NSimulations;
    error_velocity_RMS_offset_avg = sum(error_velocity_RMS_offset)/NSimulations;
    
    fprintf("The simulation had the following conditions:\n");
    if disable_canyon == 1
        fprintf("Canyon scenario disabled\n");
    else
        fprintf("Canyon scenario enabled\n");
    end
    if disable_ionosphere == 1
        fprintf("Ionospheric perturbations disabled\n");
    else
        fprintf("Ionospheric perturbations enabled\n");
    end
    
    fprintf("The noise sigma is: %f\nThe number of satellites to use for the PDOP minimization is %d\n",noise_sigma, N_PDOP);
    
    fprintf("The position RMS error for %d iterations is %f\n",NSimulations, error_position_RMS_avg);
    fprintf("The velocity RMS error for %d iterations is %f\n",NSimulations,error_velocity_RMS_avg);
    fprintf("The position RMS error for %d iterations (ignoring the first position estimation) is %f\n",NSimulations,error_position_RMS_offset_avg);
    fprintf("The velocity RMS error for %d iterations (ignoring the first velocity estimation) is %f\n",NSimulations,error_velocity_RMS_offset_avg);
end

