function [position,velocity] = update_position(time)
%UPDATE_POSITION Computes the aircraft position, velocity and acceleration
%at each time step in the N-E plane
%   This function receives time (in seconds)as input, and updates the aircraft
%   position, velocity and acceleration accordingly. The function only
%   needs to receive time since the aircraft motion is assumed perfect,
%   specified by the project handout. This motion is split into 4 segments
%   according to the simulation time. More information is present in the
%   report.

if time <=100
    %Conditions for the first leg of the trajectory
    position_initial = [0,0];
    velocity_initial = 50;
    acceleration = 1;
    direction = 90;
    
    [position,velocity]=linear_motion(time,direction,position_initial,velocity_initial,acceleration);
    
elseif time <=150
    %Conditions for the second leg of the trajectory
    position_initial = [0,10000];
    velocity_initial = 150;
    time=time-100;
    
    [position,velocity]=circular_motion(time,position_initial,velocity_initial);
    
elseif time <=200
    %Conditions for the third leg of the trajectory
    position_initial = [-1591.54943091895,8408.45056908105];
    velocity_initial = 150;
    direction = 0;
    time=time-150;
    
    [position,velocity]=linear_motion(time,direction,position_initial,velocity_initial);

    
elseif time <=250
    %Conditions for the fourth leg of the trajectory
    position_initial = [5908.45056908105,8408.45056908105];
    velocity_initial = 150;
    direction = 0;
    time=time-200;
    
    [position,velocity]=linear_motion(time,direction,position_initial,velocity_initial);
    
    
else
    fprintf("Time out of bounds\n");
end

end

