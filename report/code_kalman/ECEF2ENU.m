function [ENU] = ECEF2ENU(xyz,lat,lon)

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);

lat=deg2rad(lat);
lon=deg2rad(lon);

pi_2=pi/2;


rotx=[1,0,0;
      0,cos(pi_2-lat),-sin(pi_2-lat);
      0,sin(pi_2-lat),cos(pi_2-lat)];
  
rotz=[cos(lon+pi_2),-sin(lon+pi_2),0;
      sin(lon+pi_2),cos(lon+pi_2),0;
      0,0,1];
  
ENU=[x,y,z]*rotz*rotx;

end

