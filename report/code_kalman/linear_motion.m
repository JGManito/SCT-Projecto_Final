function [pos,v]=linear_motion(time,direction,position_initial,velocity_initial,acceleration)
%UNTITLED This function computes the position and velocity of a point in
%linear motion according to time and direction
%   This function computes the position and velocity of a point in linear 
%   motion using the formulas x=v0(t)+(a0/2)t^2 and v=a0*t, and their y 
%   counterpart. The input variable direction is the angle with regards to 
%   the x axis in degrees.
%
%   The input variable acceleration can be ommited. In this case the 
%   function computes the position as a motion at constant velocity, 
%   instead of motion at constant acceleration.
%
%   The input variable position_initial is a vector [x0 y0].
%
%   The output is always relative initial time. As such, the time variable 
%   must be offset prior to input

x0=position_initial(1);
y0=position_initial(2);

switch nargin
    case 5
        %Calculate acceleration and velocity components
        ax=acceleration*cosd(direction);
        ay=acceleration*sind(direction);
        
        vx0=velocity_initial*cosd(direction);
        vy0=velocity_initial*sind(direction);
        
        %For X
        vx=vx0 + ax*time;
        x=x0 + vx0*time + 0.5*ax*(time^2);
        
        %For Y
        vy=vy0 + ay*time;
        y=y0 + vy0*time + 0.5*ay*(time^2);
    
    case 4
        %Calculate velocity components
        vx0=velocity_initial*cosd(direction);
        vy0=velocity_initial*sind(direction);
        
        %For X
        vx=vx0;
        x=x0 + vx0*time;
        
        %For Y
        vy=vy0;
        y=y0 + vy0*time;
        
        
    otherwise
        fprintf("Error: Function 'linear_motion' doesn't have enough arguments\n");
        x=[];
        y=[];
        vx=[];
        vy=[];
end
   
pos=[x,y];
v=[vx,vy];

end
