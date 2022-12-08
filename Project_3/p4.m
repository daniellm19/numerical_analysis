clear all; close all; clc

ns = 20:10:100;
ms = 20:10:100;
reference_w = poisson(0,2,0,2,10,10,2,false);
too_long = ['m', 'n'];
too_much_error = ['m', 'n'];
good_values = ['m', 'n'];

for n = ns
    for m = ms
        tic;
        w = poisson(0,2,0,2,m,n,2,false);
        duration = toc
        %max_deviation = deviation(reference_w, w);
        if duration >= 0.5
            too_long = cat(1, too_long, [m, n]);
        end
        %if max_deviation >= 0.01
        %    too_much_error = [too_much_error; {'m': m, 'n': n, 'Max error': max_deviation}];
        %end
        %if duration < 0.5 && max_deviation < 0.01
        %    good_values = [good_values; {'m': m, 'n': n, 'Time': duration, 'Max error': max_deviation}];
        %end
    end
end

too_long
