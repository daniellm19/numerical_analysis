clear all; close all; clc;

tol = 0.001;
allowed_temp = 100;
P = []; % Max power for given K

for K = 1:0.1:5
    a = 0; % Least power
    b = 15; % Most power
    if sign(max_temp(a, allowed_temp, K))*sign(max_temp(b, allowed_temp, K)) >= 0
        error('f(a)f(b)<0 not satisfied!')
    end
    fa = max_temp(a, allowed_temp, K);
    while (b-a)/2>tol
        c=(a+b)/2;
        fc=max_temp(c, allowed_temp, K);
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
    xc=(a+b)/2;               %new midpoint is best estimate
    P = cat(1, P, xc);
end

K = 1:0.1:5;

plot(K, P)
xlabel('K [W/(cm°C)]', 'LineWidth',3)
ylabel('P [W]')