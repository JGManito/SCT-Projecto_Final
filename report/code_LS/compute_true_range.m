function [true_range] = compute_true_range(satellite_positions,receiver_position)

for i=1:size(satellite_positions,1)
    true_range(i,:)=norm(satellite_positions(i,:)-receiver_position);
end

end
