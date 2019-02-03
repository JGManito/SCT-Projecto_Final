%--------------------------------------------------------%
%           SCT - Projecto Final                         %
%--------------------------------------------------------%
%           Projecto 1                                   %
%                                                        %
%Pedro Afonso                                            %
%Joao Manito                                             %
%--------------------------------------------------------%
%Funcao main para a computa��o das posi��es usando Least %
%Squares                                                 %
%--------------------------------------------------------%

clear;
clc;

%Load ephemeris data
%Data used: Yuma week 1014 - 319488
%yuma_dataset=read_yuma("almanac.yuma.week1014.319488.txt");
yuma_dataset=read_yuma("almanac.yuma.week0961.147456.txt");

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

%Condi��es iniciais da filtragem de Kalman
x_predz=0;
vxi=0;
y_predz=0;
vyi=0;
z_predz=0;
vzi=0;
DeltaTz=0;
DeltaTponto=0;

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
disp(satellite_positions);

%Compute the true range to the receiver
true_range = compute_true_range(satellite_positions,pos_real_current_ECEF(t_sim,:));

%Determine satellites in view using the initial estimate for position
mask_angle = 10;
[satellite_positions_visible,satellite_visible,true_range_visible] = filter_satellites_visible(satellite_positions,satellite,true_range,mask_angle,position_estimate);

%Compute the pseudorange to the receiver
pseudorange = compute_pseudoranges(true_range);



%------Receiver Motion Simulation------%
z(:,1) = pseudorange_visible;
n = size(pseudorange_visible,1);

x_sat(:,1) = satellite_positions_visible(:,1);
y_sat(:,1) = satellite_positions_visible(:,2);
z_sat(:,1) = satellite_positions_visible(:,3);

%%Extended Kalman Apllied
%Ganho K da primeira itera��o

H(:,:,1)=dobs(n,x_sat,y_sat,z_sat,x_predz,y_predz,z_predz);
R(:,:,1)=varpseudo(sigma_R,n);
K(:,:,1)= gain(P_kz,H(:,:,1),R(:,:,1));


%definir h para para o estimate update
h(:,1)= obs(n,x_sat,y_sat,z_sat,x_predz,y_predz,z_predz,DeltaTz,c);


%Estimate update da primeira itera��o
states_estim(:,:,1) = update(states_predz,K(:,:,1),h(:,:,1),z(:,1));

%Error Covariance update

P_estim(:,:,1) = cov_update(K(:,:,1),H(:,:,1),P_kz,R(:,:,1),n);

%Prediction states and Covariance

states_predict(:,1)=predict_states(A,states_estim(:,:,1));

P_predict(:,:,1)= predict_cov(A,P_estim(:,:,1),Q);

position_estimate_WGS84 = XYZ2LLH([states_estim(1,1,1),states_estim(3,1,1),states_estim(5,1,1)]);

N_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,1);
E_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,2);

N_XY(t_sim)=pos_real_current_xy(t_sim,2);
E_XY(t_sim)=pos_real_current_xy(t_sim,1);

N_LLH(t_sim)=pos_real_current_WGS84(t_sim,1);
E_LLH(t_sim)=pos_real_current_WGS84(t_sim,2);


%Repetir Processo para 250 segundos
for t_sim=2:1:251 %Total trajectory time is 250s, with 1Hz sample rate
    
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
    disp(satellite_positions);
    
    %Compute the true range to the receiver
    true_range = compute_true_range(satellite_positions,pos_real_current_ECEF(t_sim,:));
    
    %Compute the pseudorange to the receiver
    pseudorange = compute_pseudoranges(true_range);
    
    
    
    %------Receiver Motion Simulation------%
    %Determine satellites in view using the initial estimate for position
    mask_angle = 10;
    [satellite_positions_visible,satellite_visible,pseudorange_visible,true_range_visible] = filter_satellites_visible(satellite_positions,satellite,pseudorange,true_range,mask_angle,position_estimate);
    
    z(:,t_sim) = pseudorange_visible;
    n = size(pseudorange_visible,1);
    x_sat(:,t_sim) = satellite_positions_visible(:,1);
    y_sat(:,t_sim) = satellite_positions_visible(:,2);
    z_sat(:,t_sim) = satellite_positions_visible(:,3);
    
    %Ganho K da primeira itera��o
    
    H(:,:,t_sim)=dobs(n,x_sat,y_sat,z_sat,states_predict(1,t_sim-1),states_predict(3,t_sim-1),states_predict(5,t_sim-1));
    R(:,:,t_sim)=varpseudo(sigma_R,n);
    K(:,:,t_sim)= gain(P_predict(:,:,t_sim-1),H(:,:,t_sim),R(:,:,t_sim));
    
    
    %Definir h para o estimate update
    
    h(:,t_sim) = obs(n,x_sat(:,t_sim),y_sat(:,t_sim),z_sat(:,t_sim),states_predict(1,t_sim-1),states_predict(3,t_sim-1),states_predict(5,t_sim-1),states_predict(6,t_sim-1),c);
    
    %Estimate update
    
    states_estim(:,:,t_sim) = update(states_predict(:,t_sim-1),K(:,:,t_sim),h(:,t_sim),z(:,t_sim));
    
    %Error Covariance update
    
    P_estim(:,:,t_sim) = cov_update(K(:,:,t_sim),H(:,:,t_sim),P_predict(:,:,t_sim-1),R(:,:,t_sim),n);
    
    
    %Prediction states and Covariance
    
    states_predict(:,t_sim)=predict_states(A,states_estim(:,:,t_sim));
    
    P_predict(:,:,t_sim)= predict_cov(A,P_estim(:,:,t_sim),Q);
    
    
    position_estimate_WGS84(t_sim,:) = XYZ2LLH([states_estim(1,1,t_sim),states_estim(3,1,t_sim),states_estim(5,1,t_sim)]);
    
    N_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,1);
    E_estimate_LLH(t_sim)=position_estimate_WGS84(t_sim,2);
    
    N_XY(t_sim)=pos_real_current_xy(t_sim,2);
    E_XY(t_sim)=pos_real_current_xy(t_sim,1);
    
    N_LLH(t_sim)=pos_real_current_WGS84(t_sim,1);
    E_LLH(t_sim)=pos_real_current_WGS84(t_sim,2);
    
    
end

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


%Matriz K - Fun��o ganho
function [K_k] = gain(P_k,H_k,R_k)

H_kT=transpose(H_k);

K_k = P_k*H_kT*inv((H_k*P_k*H_kT)+R_k);

end

%Estimate update
function x_estimn = update(states_pred,K_k,h_k,z_k)

x_estimn = states_pred + K_k*(z_k-h_k);

end

%matriz h das observa�oes usada no estimate update

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
