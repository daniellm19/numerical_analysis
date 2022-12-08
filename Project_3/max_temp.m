function w_max = max_temp(power, allowed_temp, K)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
w = poisson2(0, 4, 0, 4, 90, 40, 1, 3, power, K, 0.005, false);
w_max = max(max(w)) - allowed_temp;
end

