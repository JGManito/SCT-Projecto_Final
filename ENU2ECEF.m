function [xyz] = ENU2ECEF(ENU,lat,lon)
%ENU2ECEF Summary of this function goes here
%   Detailed explanation goes here

E=ENU(:,1);
N=ENU(:,2);
U=ENU(:,3);

lat=deg2rad(lat);
lon=deg2rad(lon);


pi_2=(pi/2)*ones(size(ENU,1),1);

rotx=[1,0,0;
      0,cos(lat-pi_2),-sin(lat-pi_2);
      0,sin(lat-pi_2),cos(lat-pi_2)];
  
rotz=[cos(-lon-pi_2),-sin(-lon-pi_2),0;
      sin(-lon-pi_2),cos(-lon-pi_2),0;
      0,0,1];
  

xyz=[E,N,U]*rotx*rotz;


end

