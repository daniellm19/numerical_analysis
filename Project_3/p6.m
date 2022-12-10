clear all; close all; clc
%w=poisson2(xl,xr,yb,yt,M,N,L_start,L_stop,P, K, H, delta, plot)

max_heat = [];
for L_start = 0:0.1:2
    w = poisson(0, 4, 0, 4, 90, 40, L_start, L_start + 2, 5, 1.68, 0.005, 0.1, false);
    max_heat = cat(1, max_heat, [L_start,max(max(w))]);
end

max_heat

L = 0:0.1:2;
plot(L, max_heat(:,2))