%% Problem 3
clear all; close all; clc;

tol = 0.01;
allowed_temp = 100;
P = []; % Max power for given H
K = 1.68;
H = 0.005;


for delta = 0.1:0.1:1
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
    xc=(a+b)/2 %new midpoint is best estimate
    P = cat(1, P, xc);
end

delta = 0.1:0.1:1;

plot(delta, P)
xlabel('\delta [cm]', 'LineWidth',3)
ylabel('P [W]')

%% Problem 4
clear all; close all; clc

p = 5;
K = 1.68;
H = 0.005;
delta = 0.1;
ns = 10:10:90;
ms = 10:10:90;
reference_w = poisson(0, 2, 0, 2, 100, 100, 0, 2, p, K, H, delta, false);
too_long = [];
too_much_error = [];
good_values = [];

for n = ns
    for m = ms
        tic;
        w = poisson(0, 2, 0, 2, m, n, 0, 2, p, K, H, delta,false);
        duration = toc;
        max_deviation = max([abs(w(1,1) - reference_w(1,1)), abs(w(end,1) - reference_w(end,1)), abs(w(1,end) - reference_w(1,end)), abs(w(end,end) - reference_w(end,end))]);
        if duration >= 0.5
            too_long = cat(1, too_long, [m, n]);
        end
        if max_deviation >= 0.01
            too_much_error = cat(1, too_much_error, [m, n]);
        end
        if duration < 0.5 && max_deviation < 0.01
            good_values = cat(1, good_values, [m,n]);
        end
    end
end

too_long
too_much_error
good_values

%% Problem 5
clear all; close all; clc

w = poisson(0, 4, 0, 4, 90, 40, 0, 2, 5, 1.68, 0.005, 0.1, true)

%% Problem 6
clear all; close all; clc

max_heat = [];
for L_start = 0:0.1:2
    w = poisson(0, 4, 0, 4, 90, 40, L_start, L_start + 2, 5, 1.68, 0.005, 0.1, false);
    max_heat = cat(1, max_heat, [L_start,max(max(w))]);
end

max_heat

L = 0:0.1:2;
plot(L, max_heat(:,2))

%% Problem 7
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

%% Problem 8

clear all; close all; clc;

tol = 0.01;
allowed_temp = 100;
P = []; % Max power for given K
H = 0.005;
delta = 0.1;

for K = 1:0.1:5
    a = 0; % Least power
    b = 15; % Most power
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
    xc=(a+b)/2;               %new midpoint is best estimate
    P = cat(1, P, xc);
end

K = 1:0.1:5;

plot(K, P)
xlabel('K [W/(cm°C)]', 'LineWidth',3)
ylabel('P [W]')

%%  Problem 9

%% Problem 10
clear all; close all; clc;
w = poissondoubleP(0,6,0,6,50,50,2,4,1,5,2,7,1.68,0.005,true)

%% Problem 11 - variable H
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

%% Problem 12 - Variable delta
clear all; close all; clc;

tol = 0.01;
allowed_temp = 100;
P = []; % Max power for given H
K = 1.68;
H = 0.005;


for delta = 0.1:0.1:1
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
    xc=(a+b)/2 %new midpoint is best estimate
    P = cat(1, P, xc);
end

delta = 0.1:0.1:1;

plot(delta, P)
xlabel('\delta [cm]', 'LineWidth',3)
ylabel('P [W]')
