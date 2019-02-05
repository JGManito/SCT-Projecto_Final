function [az,el] = ENU2AzEl(ENU)

E=ENU(:,1);
N=ENU(:,2);
U=ENU(:,3);

az=rad2deg(atan2(E,N));

el=rad2deg(atan2(U,(sqrt(N.^2+E.^2))));

end

