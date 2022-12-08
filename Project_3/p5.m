clear all; close all; clc
%w=poisson2(xl,xr,yb,yt,M,N,L_start,L_stop,P, K, H, plot)

w = poisson2(0, 4, 0, 4, 90, 40, 0, 2, 5, 1.68, 0.005, true)