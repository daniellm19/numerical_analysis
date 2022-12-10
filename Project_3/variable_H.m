clear all; close all; clc;

tol = 0.01;
allowed_temp = 100;
P = []; % Max power for given H
K = 1.68;
delta = 0.01;

for H = 0.005:0.1:3
    a = 0; % Least power
    b = 300; % Most power
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
    xc=(a+b)/2; %new midpoint is best estimate
    P = cat(1, P, xc);
end

H = 0.005:0.1:3;

plot(H, P)
xlabel('H [W/(cm^2°C)]', 'LineWidth',3)
ylabel('P [W]')