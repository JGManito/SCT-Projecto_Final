function [llh] = flat2llh(position_xy, position_prev_xy, position_prev_llh)

%WGS84 variables
a = 6378137;
f = 1/298.257223563;
	
angle=0;

NE = [cosd(angle) -sind(angle); sind(angle) cosd(angle)]*[position_xy(2);position_xy(1)];

RN = a/sqrt(1-(2*f - (f^2))*sind(position_prev_llh(1))^2);

RM = RN*((1-(2*f-f^2))/(1-(2*f-f^2)*sind(position_prev_llh(1)^2)));

dN=NE(1)-position_prev_xy(2);
dE=NE(2)-position_prev_xy(1);

dlat=atan2d(1,RM)*dN;
lat=position_prev_llh(1)+dlat;

dlon=atan2d(1,(RN*cosd(lat)))*dE;
lon=position_prev_llh(2)+dlon;

%From the trajectory definition
h=2000;

llh=[lat,lon,h];
end

