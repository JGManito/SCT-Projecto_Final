function [G] = compute_G(satellite_positions,initial_estimate)
%COMPUTE_G Summary of this function goes here
%   Detailed explanation goes here

NSats = size(satellite_positions,1);

G=[];

for i=1:NSats
    [XYZ] = satellite_positions(i,:);
    [xyz] = initial_estimate(:);
    
    X = XYZ(1);
    Y = XYZ(2);
    Z = XYZ(3);
    
    x = xyz(1);
    y = xyz(2);
    z = xyz(3);
    
    D = norm(satellite_positions(i,:) - initial_estimate(:));
    
    dpdx = -(X-x)/D;
    dpdy = -(Y-y)/D;
    dpdz = -(Z-z)/D;
    
    G=[G;dpdx  dpdy  dpdz  1];

%     dpdx = (X-x)/D;
%     dpdy = (Y-y)/D;
%     dpdz = (Z-z)/D;
%     
%     G=[G;dpdx  dpdy  dpdz  -1];

end
    


end

