function [satellite_positions_visible,satellite_visible,pseudorange_visible,true_range_visible] = filter_satellites_visible(satellite_positions,satellite,pseudorange,true_range,mask_angle,position_estimate)
%FILTER_SATELLITES_VISIBLE Summary of this function goes here
%   Detailed explanation goes here

NSats = size(satellite,2);
j=1;

for i=1:NSats
    %Convert satellite position to ENU, relative to the position estimate
    position_ENU = ECEF2ENU(satellite_positions(i,:),position_estimate(1),position_estimate(2));
    
    %Obtain Elevation angle
    [~,el] = ENU2AzEl(position_ENU);
    satellite(i).elevation=el;
    
    if el >= mask_angle
        satellite_positions_visible(j,:) = satellite_positions(i,:);
        satellite_visible(j) = satellite(i);
        pseudorange_visible(j,:) = pseudorange(i,:);
        true_range_visible(j,:) = true_range(i,:);
        j=j+1;
    end
        




end

