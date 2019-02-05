function [LLH] = XYZ2LLH(xyz)
%XYZ2LLH converts ECEF coordinates to LLH
%   This function converts XYZ coordinates in the ECEF frame to Latitude,
%   Longitude and Altitude coordinates in the WGS84 coordinate system. The
%   formulas used where obtained from Kaplan, 2018 (chapter 2.2.5)
x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);

a=6378137;
f=1/298.257223563;
b=a*(1-f);
e=sqrt(1-(b^2 / a^2));
e_dash=(a/b)*e;


if x>=0
    lon=atand(y/x);
elseif (x<0 && y>=0)
    lon=180+atand(y/x);
elseif (x<0 && y<0)
    lon=-180+atand(y/x);
end

p=sqrt(sum([x,y].^2,2));
u=atan2d(z*a,p*b);

prev_iter=1;

%iteration loop
while (abs(tand(u)-prev_iter) > eps)
    cos_u = sqrt(1/(1+(tand(u)^2)));
    sin_u = sqrt(1-cos_u^2);
    lat=atan2d(z+(e_dash^2)*b*(sin_u^3),p-(e^2)*a*(cos_u^3));
    u=atan2d(b*tand(lat),a);
    prev_iter = tand(u);
end

N = a/sqrt(1-(e^2)*sind(lat)^2);

h=(p/cosd(lat))-N;

LLH=[lat,lon,h];
        
 
end

