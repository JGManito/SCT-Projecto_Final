function [satellite_positions_visible,satellite_visible,true_range_visible,elevation,NView] = filter_satellites_visible_LS(satellite_positions,satellite,true_range,mask_angle,position_estimate,flag)

NSats = size(satellite,2);
j=1;

LLH_initial_estimate = XYZ2LLH(position_estimate);

for i=1:NSats
    %Convert satellite position to ENU, relative to the position estimate
    position_ENU = ECEF2ENU(satellite_positions(i,:),LLH_initial_estimate(1),LLH_initial_estimate(2));
    
    %Obtain Elevation angle
    [~,el] = ENU2AzEl(position_ENU);
    satellite(i).elevation=el;
    
    if el >= mask_angle
        satellite_positions_visible(j,:) = satellite_positions(i,:);
        satellite_visible(j) = satellite(i);
        true_range_visible(j,:) = true_range(i,:);
        elevation(j) = el;
        j=j+1;
    end
    canyon_sats_index = [];
end

NView = j-1;

if flag == 1
    %Save the filtered results to temporary variables
    positions_temp(:,:) = satellite_positions_visible(:,:);
    satellite_visible_temp(:) = satellite_visible(:);
    true_range_visible_temp(:,:) = true_range_visible(:,:);
    elevation_temp(:) = elevation(:);
    
    %Find the 4 highest elevations
    [~,canyon_sats_index] = maxk(elevation_temp,NSats);
    
    %Reset the output variables
    satellite_positions_visible = [];
    satellite_visible (NSats+1:end)= [];
    true_range_visible= [];
    elevation = [];
    
    %Output only the satellites with the 4 highest elevations
    for j=1:4
        [~,i]=min(canyon_sats_index);
        satellite_positions_visible(j,:) = positions_temp(canyon_sats_index(i),:);
        satellite_visible(j) = satellite_visible_temp(canyon_sats_index(i));
        true_range_visible (j,1)= true_range_visible_temp(canyon_sats_index(i));
        elevation(j) = elevation_temp(canyon_sats_index(i));
        
        canyon_sats_index(i)=[];
    end
end

end


