function [satellite_positions_visible,satellite_visible,true_range_visible,elevation] = filter_satellites_visible(satellite_positions,satellite,true_range,mask_angle,position_estimate,flag)
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
        %pseudorange_visible(j,:) = pseudorange(i,:);
        true_range_visible(j,:) = true_range(i,:);
        elevation(j) = el;
        j=j+1;
    end
    
end

if flag == 1
    %Save the filtered results to temporary variables
    positions_temp(:,:) = satellite_positions_visible(:,:);
    satellite_visible_temp(:) = satellite_visible(:);
    %pseudorange_visible_temp(:,:) = pseudorange_visible(:,:);
    true_range_visible_temp(:,:) = true_range_visible(:,:);
    elevation_temp(:) = elevation(:);
    
    %Find the 4 highest elevations
    canyon_sats_index = maxk(elevation_temp,4);
    
    %Reset the output variables
    satellite_positions_visible = [];
    satellite_visible= [];
    %pseudorange_visible= [];
    true_range_visible= [];
    elevation = [];
    
    %Output only the satellites with the 4 highest elevations
    satellite_positions_visible(:,:) = [positions_temp(canyon_sats_index(1),:);positions_temp(canyon_sats_index(2),:);positions_temp(canyon_sats_index(3),:);positions_temp(canyon_sats_index(4),:)];
    satellite_visible(:) = [satellite_visible_temp(canyon_sats_index(1)),satellite_visible_temp(canyon_sats_index(2)),satellite_visible_temp(canyon_sats_index(3)),satellite_visible_temp(canyon_sats_index(4))];
    %pseudorange_visible(:,:) = [pseudorange_visible_temp(canyon_sats_index(1),:),pseudorange_visible_temp(canyon_sats_index(2),:),pseudorange_visible_temp(canyon_sats_index(3),:),pseudorange_visible_temp(canyon_sats_index(4),:)];
    true_range_visible (:,:)= [true_range_visible_temp(canyon_sats_index(1),:),true_range_visible_temp(canyon_sats_index(2),:),true_range_visible_temp(canyon_sats_index(3),:),true_range_visible_temp(canyon_sats_index(4),:)];
    elevation(:) = [elevation_temp(canyon_sats_index(1)),elevation_temp(canyon_sats_index(2)),elevation_temp(canyon_sats_index(3)),elevation_temp(canyon_sats_index(4))];
end

end


