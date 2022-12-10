function w_max = max_temp(power, allowed_temp, K, H, delta)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
w = poisson(0, 4, 0, 4, 90, 40, 1, 3, power, K, H, delta, false);
w_max = max(max(w)) - allowed_temp;
end

