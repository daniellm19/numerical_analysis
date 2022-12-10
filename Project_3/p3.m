clear all; close all; clc

p = 5;
K = 1.68;
H = 0.005;
delta = 0.1;

%w=poisson(xl,xr,yb,yt,M,N,L_start,L_stop,p, K, H, delta, plot)

w = poisson(0, 2, 0, 2, 10, 10, 0, 2, p, K, H, delta, true)


