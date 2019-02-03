function [satellite_parameters] = compute_orbital_parameters(yuma_data)
%COMPUTE Summary of this function goes here
%   Detailed explanation goes here

NSats=size(yuma_data,2);


for i=1:NSats
%Reasign parsed data into local variables for easier reading
ID = yuma_data(i).id;
e = yuma_data(i).eccentricity;
toe = yuma_data(i).time;
inc = yuma_data(i).inclination;
OmegaDot = yuma_data(i).omega_dot;
sqrtA = yuma_data(i).sqrt_A;
Omega0 = yuma_data(i).omega_zero;
argP = yuma_data(i).arg_perigee;
M0 = yuma_data(i).mean_anomaly;
af0 = yuma_data(i).af0;
af1 = yuma_data(i).af1;
WN = yuma_data(i).week;

%Change almanac data to other forms, more suitable for the present project
toe = (WN + 1024) * 7 * 24 * 3600 + toe; %Obtain toe since the start of GPS time
inc = rad2deg(inc); %Convert inclination to degrees
OmegaDot = rad2deg(OmegaDot); %Convert rate of change of the ascending node to degrees
A = sqrtA^2; %Use semi-major axis instead of it's square root
Omega0 = rad2deg(Omega0); %Convert longitude of the ascending node at the begining of the GPS week to degrees
argP = rad2deg(argP); %Convert the argument of perigee to degrees;
M0 = rad2deg(M0); %Convert the mean anomaly to degrees


%Compute variables needed to define the satellite orbital position
R=26559800;     %Orbital radius
inc=55;           %Orbit inclination angle
T=43082;        %Orbital period

%Assign all the variables to output structure
satellite_parameters(i).ID = ID;
satellite_parameters(i).R = R;
satellite_parameters(i).T = T;
satellite_parameters(i).e = e;
satellite_parameters(i).toe = toe;
satellite_parameters(i).inc = inc;
satellite_parameters(i).OmegaDot = OmegaDot;
satellite_parameters(i).A = A;
satellite_parameters(i).Omega0 = Omega0;
satellite_parameters(i).argP = argP;
satellite_parameters(i).M0 = M0;
satellite_parameters(i).af0 = af0;
satellite_parameters(i).af1 = af1;
satellite_parameters(i).WN = WN;
end

end

