function [positions,satellite] = compute_satellite_position(satellite,time)

NSats = size(satellite,2);

for i=1:NSats
    %Obtain simulation time, referenced to the start of GPS time
    WN = satellite(i).WN;
    time = 7 * 24 * 3600 + 432000 + time; %Set reference time to Feb 1 2019 00:00 UTC
    
    dt = time - satellite(i).toe;
    if dt > 302400
        dt = dt-604800;
    end
    if dt < -302400
        dt = dt + 604800;
    end
    
    
    %Compute argument of latitude and longitude of ascending node
    argLat = (satellite(i).argP) + (360/(satellite(i).T))*dt;
    Omega = satellite(i).Omega0 - (180/satellite(i).T)*time;
    
    
    %Compute coordinates in ECEF
    R = satellite(i).R;
    alpha = satellite(i).inc;
    XYZ = [(R*cosd(argLat)*cosd(Omega) - R*sind(argLat)*sind(Omega)*cosd(alpha));...
        R*(cosd(argLat)*sind(Omega) + sind(argLat)*cosd(Omega)*cosd(alpha));...
        R*(sind(argLat)*sind(alpha))];
    
    %Save positions to output variable
    positions(i,:) = XYZ';
    
    %Add new orbital parameters to satellite variable
    satellite(i).argLat = argLat;
    satellite(i).Omega = Omega; 
end

end

