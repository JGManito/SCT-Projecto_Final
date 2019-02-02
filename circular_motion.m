function [pos,v]=circular_motion(time,position_initial,velocity_initial)
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


T=(4/3)*50;
omega=(2*pi)/T;
R=150/omega;

x_diff = position_initial(1) - R;
y_diff = position_initial(2) - 0;

%theta_initial=atan2d(y_,x_);

theta=omega*time;
x_ = R*cos(theta);
y_ = R*sin(theta);


x=x_diff +x_;
y=y_diff +y_;

vx=velocity_initial*cos(theta+0.5*pi);
vy=velocity_initial*sin(theta+0.5*pi);

pos=[x,y];
v=[vx,vy];



end