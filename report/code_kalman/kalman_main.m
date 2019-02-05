%--------------------------------------------------------%
%           SCT - Projecto Final                         %
%--------------------------------------------------------%
%           Projecto 1                                   %
%                                                        %
%Pedro Afonso                                            %
%Joao Manito                                             %
%--------------------------------------------------------%
%Funcao main para a computaï¿½ï¿½o das posiï¿½ï¿½es usando Kalman%                                                 %
%--------------------------------------------------------%

clear;
clc;

%Define problem variables
mask_angle = 10;
noise_sigma = sqrt(10);
N_PDOP = 5;
disable_canyon = 0;
disable_noise = 0;
disable_ionosphere = 0;
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
    
    %Definir constantes
    c=3*(10^8);
    deltat= 1;
    h_mdois=2*(10^-20);
    h_zero=2*(10^-19);
    sigma_R=7.11;
    
    %Condiï¿½ï¿½es iniciais da filtragem de Kalman
    x_predz=0;
    vxi=0;
    y_predz=0;
    vyi=0;
    z_predz=0;
    vzi=0;
    DeltaTz=0;
    DeltaTponto=0;
    states_estim = [];
    states_predict = [];
    P_predict = [];
    
    states_predz=[x_predz,vxi,y_predz,vyi,z_predz,vzi,DeltaTz,DeltaTponto]';
    
    q_v = 10;
    
    %Allan Variance parameters
    q_f = 2*(pi^2)*h_mdois;
    q_phi = h_zero/2;
    
    
    a=(q_phi*deltat + (q_f*(deltat^3)/3))*(c^2);
    b=(q_f*(deltat^2)*(c^2)/2);
    c= q_f*(c^2)*deltat;
    
    
    
    A = [1 deltat 0 0 0 0 0 0 ;
        0 1 0 0 0 0 0 0 ;
        0 0 1 deltat 0 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 0 0 deltat 1 0 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 deltat 1;
        0 0 0 0 0 0 0 1];
    
    
    Q = [q_v*(deltat^3)/3 q_v*(deltat^2)/2 0 0 0 0 0 0;
        q_v*(deltat^2)/2 q_v*deltat 0 0 0 0 0 0;
        0 0 q_v*(deltat^3)/3 q_v*(deltat^2)/2 0 0 0 0;
        0 0 q_v*(deltat^2)/2 q_v*deltat 0 0 0 0;
        0 0 0 0 q_v*(deltat^3)/3 q_v*(deltat^2)/2 0 0;
        0 0 0 0 q_v*(deltat^2)/2 q_v*deltat 0 0;
        0 0 0 0 0 0 a b ;
        0 0 0 0 0 0 b c];
    
    %Covariance Matrix Errors initial
    
    P_kz=[10^5 0 0 0 0 0 0 0;
        0 10^5 0 0 0 0 0 0
        0 0 10^5 0 0 0 0 0
        0 0 0 10^5 0 0 0 0
        0 0 0 0 10^5 0 0 0
        0 0 0 0 0 10^5 0 0
        0 0 0 0 0 0 10^5 0
        0 0 0 0 0 0 0 10^5];
    
    
    %Start the simulation
    t_sim = 1;
    fprintf("Iteration #%d, time =%d\n",iteration,t_sim-1);
    
    %------Aircraft Motion Simulation------%
    pos_real_prev_xy = pos_real_current_xy(1,:);
    pos_real_prev_WGS84 = pos_real_current_WGS84(1,:);
    pos_real_prev_ECEF = LLH2XYZ(pos_real_current_WGS84(1,:));
    
    %Update position in the XY plane
    [pos_real_current_xy(t_sim,:),velocity_real_current_xy(t_sim,:)]=update_position(t_sim-1);
    
    %Compute the new position in WGS84 and ECEF frames
    pos_real_current_WGS84(t_sim,:) = flat2llh(pos_real_current_xy(t_sim,:) , pos_real_prev_xy , pos_real_prev_WGS84);
    pos_real_current_ECEF(t_sim,:) = LLH2XYZ(pos_real_current_WGS84(t_sim,:));
    
    
    
    
    %------Satellite Motion Simulation------%
    %Get satellite positions in ECEF
    satellite = compute_orbital_parameters(yuma_dataset);
    [satellite_positions,satellite] = compute_satellite_position(satellite,t_sim-1);
    
    %Compute the true range to the receiver
    true_range = compute_true_range(satellite_positions,pos_real_current_ECEF(t_sim,:));
    
    %Determine satellites in view using the initial estimate for position
    [satellite_positions_visible,satellite_visible,true_range_visible,elevation] = filter_satellites_visible_KALMAN(satellite_positions,satellite,true_range,mask_angle,position_estimate,0);
    
    %Compute the pseudorange to the receiver
    pseudorange_visible = compute_pseudoranges(true_range_visible,elevation,noise_sigma,disable_ionosphere,disable_noise);
    
    [pseudorange_visible,satellite_positions_visible] = minimize_PDOP(satellite_positions_visible,pseudorange_visible,position_estimate,N_PDOP);
    
    
    %------Receiver Motion Simulation------%
    z = pseudorange_visible;
    n = size(pseudorange_visible,1);
    
    x_sat = satellite_positions_visible(:,1);
    y_sat = satellite_positions_visible(:,2);
    z_sat = satellite_positions_visible(:,3);
    
    %%Extended Kalman Apllied
    %Ganho K da primeira iteraï¿½ï¿½o
    
    H=dobs(n,x_sat,y_sat,z_sat,x_predz,y_predz,z_predz);
    R=varpseudo(sigma_R,n);
    K= gain(P_kz,H(:,:),R(:,:));
    
    
    %definir h para para o estimate update
    h= obs(n,x_sat,y_sat,z_sat,x_predz,y_predz,z_predz,DeltaTz,c);
    
    
    %Estimate update da primeira iteraï¿½ï¿½o
    states_estim(:,1) = update(states_predz,K(:,:),h(:),z(:));
    
    %Error Covariance update
    
    P_estim = cov_update(K(:,:),H(:,:),P_kz,R(:,:),n);
    
    %Prediction states and Covariance
    
    states_predict(:,1)=predict_states(A,states_estim(:,:,1));
    
    P_predict(:,:,1)= predict_cov(A,P_estim(:,:,1),Q);
    
    position_estimate_WGS84(1,:) = XYZ2LLH([states_estim(1,1,1),states_estim(3,1,1),states_estim(5,1,1)]);
    position_estimate_ECEF(1,:) = [states_estim(1,1,1),states_estim(3,1,1),states_estim(5,1,1)];
    
    N_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,1);
    E_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,2);
    
    N_XY(t_sim)=pos_real_current_xy(t_sim,2);
    E_XY(t_sim)=pos_real_current_xy(t_sim,1);
    
    N_LLH(t_sim)=pos_real_current_WGS84(t_sim,1);
    E_LLH(t_sim)=pos_real_current_WGS84(t_sim,2);
    
    
    %Repetir Processo para 250 segundos
    for t_sim=2:1:251 %Total trajectory time is 250s, with 1Hz sample rate
        fprintf("Iteration #%d, time =%d\n",iteration,t_sim-1);
        
        %Change to N=4 for canyon scenario
        if t_sim >= 202 && disable_canyon == 0
            flag = 1;
        else
            flag = 0;
        end
        
        %------Aircraft Motion Simulation------%
        pos_real_prev_xy = pos_real_current_xy(t_sim-1,:);
        pos_real_prev_WGS84 = pos_real_current_WGS84(t_sim-1,:);
        pos_real_prev_ECEF = LLH2XYZ(pos_real_current_WGS84(t_sim-1,:));
        
        %Update position in the XY plane
        [pos_real_current_xy(t_sim,:),velocity_real_current_xy(t_sim,:)]=update_position(t_sim-1);
        
        %Compute the new position in WGS84 and ECEF frames
        pos_real_current_WGS84(t_sim,:) = flat2llh(pos_real_current_xy(t_sim,:) , pos_real_prev_xy , pos_real_prev_WGS84);
        pos_real_current_ECEF(t_sim,:) = LLH2XYZ(pos_real_current_WGS84(t_sim,:));
        
        
        
        
        %------Satellite Motion Simulation------%
        %Get satellite positions in ECEF
        satellite = compute_orbital_parameters(yuma_dataset);
        [satellite_positions,satellite] = compute_satellite_position(satellite,t_sim-1);
        
        %Compute the true range to the receiver
        true_range = compute_true_range(satellite_positions,pos_real_current_ECEF(t_sim,:));
        
        
        %Determine satellites in view using the initial estimate for position
        [satellite_positions_visible,satellite_visible,true_range_visible,elevation,canyon_sats_index] = filter_satellites_visible_KALMAN(satellite_positions,satellite,true_range,mask_angle,position_estimate,flag);
        
        %Compute the pseudorange to the receiver
        pseudorange_visible = compute_pseudoranges(true_range_visible,elevation,noise_sigma,disable_ionosphere,disable_noise);
        
        position_estimate = LLH2XYZ(position_estimate_WGS84(t_sim-1,:));
        
        if flag == 0
            [pseudorange_visible,satellite_positions_visible] = minimize_PDOP(satellite_positions_visible,pseudorange_visible,position_estimate,N_PDOP);
        end
        
        
        
        %------Receiver Motion Simulation------%
        if flag == 1
            z = pseudorange_visible;
            n = size(pseudorange_visible,1);
            x_sat = satellite_positions_visible(:,1);
            y_sat = satellite_positions_visible(:,2);
            z_sat = satellite_positions_visible(:,3);
        else
            z = pseudorange_visible;
            n = size(pseudorange_visible,1);
            x_sat = satellite_positions_visible(:,1);
            y_sat = satellite_positions_visible(:,2);
            z_sat = satellite_positions_visible(:,3);
        end
        
        %Ganho K da primeira iteraï¿½ï¿½o
        
        H = dobs(n,x_sat,y_sat,z_sat,states_predict(1,t_sim-1),states_predict(3,t_sim-1),states_predict(5,t_sim-1));
        R = varpseudo(sigma_R,n);
        K = gain(P_predict(:,:,t_sim-1),H(:,:),R(:,:));
        
        
        %Definir h para o estimate update
        
        h = obs(n,x_sat(:),y_sat(:),z_sat(:),states_predict(1,t_sim-1),states_predict(3,t_sim-1),states_predict(5,t_sim-1),states_predict(6,t_sim-1),c);
        
        %Estimate update
        
        states_estim(:,t_sim) = update(states_predict(:,t_sim-1),K(:,:),h(:),z(:));
        
        %Error Covariance update
        
        P_estim = cov_update(K(:,:),H(:,:),P_predict(:,:,t_sim-1),R(:,:),n);
        
        
        %Prediction states and Covariance
        
        states_predict(:,t_sim)=predict_states(A,states_estim(:,t_sim));
        
        P_predict(:,:,t_sim)= predict_cov(A,P_estim(:,:),Q);
        
        position_estimate_ECEF(t_sim,:) = [states_estim(1,t_sim),states_estim(3,t_sim),states_estim(5,t_sim)];
        position_estimate_WGS84(t_sim,:) = XYZ2LLH([states_estim(1,t_sim),states_estim(3,t_sim),states_estim(5,t_sim)]);
        
        N_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,1);
        E_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,2);
        
        VX(t_sim) = states_estim(2,t_sim);
        VY(t_sim) = states_estim(4,t_sim);
        VZ(t_sim) = states_estim(6,t_sim);
        
        VABS(t_sim) = sqrt(VX(t_sim)^2 + VY(t_sim)^2 + VZ(t_sim)^2);
        
        N_XY(t_sim)=pos_real_current_xy(t_sim,2);
        E_XY(t_sim)=pos_real_current_xy(t_sim,1);
        
        N_LLH(t_sim)=pos_real_current_WGS84(t_sim,1);
        E_LLH(t_sim)=pos_real_current_WGS84(t_sim,2);
        
        
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
    
    for i=1:251
        if abs(error_position(i)) < 100 %Ignore position error until convergence
            n_conv = i;
            break;
        end
    end
    error_position_RMS_offset(iteration) = sqrt((sum(error_position(n_conv:end).^2))/(251-n_conv));
    error_velocity_RMS_offset(iteration) = sqrt((sum(error_velocity(n_conv:end).^2))/(251-n_conv));
    
    
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
    
    fprintf("The noise sigma is: %f\nThe number of satellites to use for the PDOP minimization is %d",noise_sigma, N_PDOP);
    
    fprintf("The position RMS error for %d iterations is %f\n",NSimulations, error_position_RMS_avg);
    fprintf("The velocity RMS error for %d iterations is %f\n",NSimulations,error_velocity_RMS_avg);
    fprintf("The position RMS error for %d iterations (ignoring the first position estimation) is %f\n",NSimulations,error_position_RMS_offset_avg);
    fprintf("The velocity RMS error for %d iterations (ignoring the first velocity estimation) is %f\n",NSimulations,error_velocity_RMS_offset_avg);
end


%Matriz H
function [H_k] = dobs(n,x_sat,y_sat,z_sat,x_pred,y_pred,z_pred)
for i=1:n
    
    d=sqrt((x_sat(i)-x_pred)^2+(y_sat(i)-y_pred)^2+(z_sat(i)-z_pred)^2);
    
    H_k(i,1)= -(x_sat(i)-x_pred)/d;
    H_k(i,2)= 0;
    H_k(i,3)= -(y_sat(i)-y_pred)/d;
    H_k(i,4) = 0;
    H_k(i,5)= -(z_sat(i)-z_pred)/d;
    H_k(i,6)= 0;
    H_k(i,7)= 1;
    H_k(i,8)= 0;
end
end

%Matriz R
function [R] = varpseudo(sigma_R,n)

for i=1:n
    Rvar(i)=sigma_R^2;
    
end

R=diag(Rvar);

end


%Matriz K - Funï¿½ï¿½o ganho
function [K_k] = gain(P_k,H_k,R_k)

H_kT=transpose(H_k);

K_k = P_k*H_kT*inv((H_k*P_k*H_kT)+R_k);

end

%Estimate update
function x_estimn = update(states_pred,K_k,h_k,z_k)

x_estimn = states_pred + K_k*(z_k-h_k);

end

%matriz h das observaï¿½oes usada no estimate update

function [h] = obs(n,x_sat,y_sat,z_sat,x_pred,y_pred,z_pred,DeltaT,c)

for i=1:n
    h(i,1)= sqrt((x_sat(i)-x_pred)^2+(y_sat(i)-y_pred)^2+(z_sat(i)-z_pred)^2)- c*DeltaT;
    
end
end

%Error covariance update

function [P_estim] = cov_update(K_k,H_k,P_a,R_k,n)

I = eye(8);
X = I-(K_k*H_k);
X_t = transpose(X);
K_t= transpose(K_k);
P_estim = (I-(K_k*H_k))* P_a * X_t + (K_k*R_k*K_t);

end

%Prediction states

function [x_n] = predict_states(A,states_estim)

x_n = A*states_estim;

end


%Prediction Covariance
function [P_k] = predict_cov(A,P_estim,Q)

B = transpose(A);

P_k = (A * P_estim * B) + Q;
end
