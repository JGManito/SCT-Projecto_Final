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
clf;

%Load ephemeris data
%Data used: Yuma week 1014 - 319488
yuma_dataset=read_yuma("almanac.yuma.week1014.319488.txt");

%Initialize position
pos_real_current_WGS84=[40,-9,2000];
pos_real_current_ECEF=LLH2XYZ(pos_real_current_WGS84);

%Initialize position and velocity in XY plane
pos_real_current_xy=[0,0];
velocity_real_current_xy=[0,0];

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
    
    
    
    
    
    
    N_XY(t_sim+1)=pos_real_current_xy(t_sim+1,2);
    E_XY(t_sim+1)=pos_real_current_xy(t_sim+1,1);
    
    N_LLH(t_sim+1)=pos_real_current_WGS84(t_sim+1,1);
    E_LLH(t_sim+1)=pos_real_current_WGS84(t_sim+1,2);
    
    %close all
    
    %plot(E_XY,N_XY);
    %grid on
    %daspect([1 1 1])
    
    %figure
%     plot(E_LLH,N_LLH);
%     grid on
%     daspect([1 1 1])
    
    
end

plot(E_XY,N_XY);
grid on
daspect([1 1 1])

figure
plot(E_LLH,N_LLH);
grid on
daspect([1 1 1])