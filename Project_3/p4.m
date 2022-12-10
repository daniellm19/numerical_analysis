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
