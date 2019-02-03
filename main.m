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

%Start the simulation
for t_sim=0:1:250 %Total trajectory time is 250s, with 1Hz sample rate
    
    %------Aircraft Motion Simulation------%
    %Save previous position
    if t_sim == 0
        pos_real_prev_xy = pos_real_current_xy(1,:);
        pos_real_prev_WGS84 = pos_real_current_WGS84(1,:);
        pos_real_prev_ECEF = LLH2XYZ(pos_real_current_WGS84(1,:));
    else
        pos_real_prev_xy = pos_real_current_xy(t_sim,:);
        pos_real_prev_WGS84 = pos_real_current_WGS84(t_sim,:);
        pos_real_prev_ECEF = LLH2XYZ(pos_real_current_WGS84(t_sim,:));
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
    
    %Compute the pseudorange to the receiver
    pseudorange = compute_pseudoranges(true_range);
    
    
    
    %------Receiver Motion Simulation------%
    %Determine satellites in view using the initial estimate for position
    mask_angle = 10;
    [satellite_positions_visible,satellite_visible,pseudorange_visible,true_range_visible] = filter_satellites_visible(satellite_positions,satellite,pseudorange,true_range,mask_angle,position_estimate);
    
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
    

    
end

plot(E_XY,N_XY);
grid on
daspect([1 1 1])

figure
plot(E_LLH,N_LLH,E_estimate_LLH,N_estimate_LLH);
grid on
daspect([1 1 1])