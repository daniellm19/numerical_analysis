clear all; close all; clc;
a = 5;
b = 10;
tol = 0.001;
allowed_temp = 100;
K = 1.68;
H = 0.005;
delta = 0.1;


if sign(max_temp(a, allowed_temp, K, H, delta))*sign(max_temp(b, allowed_temp, K, H, delta)) >= 0
    error('f(a)f(b)<0 not satisfied!')
end
fa = max_temp(a, allowed_temp, K, H, delta);
while (b-a)/2>tol
    c=(a+b)/2;
    fc=max_temp(c, allowed_temp, K, H, delta);
    if fc == 0              %c is a solution, done
    break
    end
    if sign(fc)*sign(fa) < 0 %a and c make the new interval
        b = c;
    else                    %c and b make the new interval
    a = c;
    fa = fc;
    end
end
xc=(a+b)/2               %new midpoint is best estimate