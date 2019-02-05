function [pseudorange] = compute_pseudoranges(true_range,elevation,noise_sigma,disable_ionosphere,disable_noise)

%Compute ionospheric influence
if disable_ionosphere == 0
    ionospheric_disturbance = 10./sind(elevation);
end

%Compute noise
for i=1:size(true_range,1)
    if disable_noise == 1
        noise(i) = 0;
    else
        noise(i) = (rand()-0.5)*2*noise_sigma;
    end
end

for i=1:size(true_range,1)
    if disable_ionosphere == 0
        pseudorange(i,1) = true_range(i) + ionospheric_disturbance(i) + noise(i);
    else
        pseudorange(i,1) = true_range(i) + noise(i);
    end
end

end

