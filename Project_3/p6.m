clear all; close all; clc
%w=poisson2(xl,xr,yb,yt,M,N,L_start,L_stop,P, K, H, plot)

max_heat = [];
for L_start = 0:0.1:2
    w = poisson2(0, 4, 0, 4, 90, 40, L_start, L_start + 2, 5, 1.68, 0.005, false);
    max_heat = cat(1, max_heat, [L_start,max(max(w))]);
end

max_heat

min(max_heat(:,2))