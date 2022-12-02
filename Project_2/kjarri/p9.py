import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos, pi, std
import matplotlib.animation as animation
import json

#Global constants
g = 9.81

def ydot(t:float, y: list, L: float, m:float):
    theta_1pp = ((m*L*pow(y[1],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m*g*sin(y[2])*cos(y[2] - y[0])) + m * L * pow(y[3],2) * sin(y[2] - y[0]) - (m+m) * g*sin(y[0]))/((m+m)*L - m*L*pow(cos(y[2] - y[0]),2))
    theta_2pp = (-1*(m*L*pow(y[3],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m+m)*(g*sin(y[0])*cos(y[2] - y[0])-L*pow(y[1], 2) * sin(y[2] - y[0])-g*sin(y[2])))/((m+m)*L - m*L*pow(cos(y[2] - y[0]), 2))
    return np.array([y[1], theta_1pp, y[3], theta_2pp])


def runge_kutta(x, n, T, L: float, m):
    h = T/n
    t = 0
    for _ in range(n):
        k1 = ydot(t, x, L, m)
        k2 = ydot(t + h/2, x + h/2 * k1, L, m)
        k3 = ydot(t + h/2, x + h/2 * k2, L, m)
        k4 = ydot(t + h, x + h * k3, L, m)
        y = np.array(x + h * (k1/6 + k2/3 + k3/3 + k4/6))
        x = y
        t += h


    return y, t, h

def random_theta():
    '''Reads the random angles from random_ini_values.json'''
    with open('random_ini_values.json') as f:
        all_ini_vals = json.load(f)
    return all_ini_vals
    
def get_error(ini_val, good_n, ns, T, L, m):
    y_0 = np.array(ini_val)
    y_good_est, t, h = runge_kutta(y_0, good_n, T, L, m)
    errors_norm = []
    for n in ns:
        y_bad_est, t, h = runge_kutta(y_0, n, T, L, m)
        errors_norm.append(np.linalg.norm(np.array([abs(y_good_est[0] - y_bad_est[0]), abs(y_good_est[1] - y_bad_est[1]), abs(y_good_est[2] - y_bad_est[2]), abs(y_good_est[3] - y_bad_est[3])])))
    slope = np.polyfit(np.log(ns), np.log(errors_norm), 1)[0]
    return errors_norm, slope




def main():
    T = 20
    good_n = 35000
    ns = np.array([100,200,400,800,1600,3200,6400])
    L = 2
    m = 1
    all_ini_vals = random_theta()
    all_errors_norm = []
    all_slopes = []
    for ini_val in all_ini_vals:
        error_norm, slope = get_error(ini_val, good_n, ns, T, L, m)
        all_errors_norm.append(error_norm)
        all_slopes.append(slope)
    print('Average:', np.average(all_slopes))
    print('Min:', np.min(all_slopes))
    print('Max:', np.max(all_slopes))
    print('Median:', np.median(all_slopes))
    print('std:', std(all_slopes))
    plt.subplot(121)
    for e in all_errors_norm:
        plt.loglog(ns, e)
    plt.ylabel('Error')
    plt.xlabel('n-value')
    plt.subplot(122)
    plt.hist(all_slopes)
    plt.ylabel('Error')
    plt.xlabel('n-value')
    plt.show()



main()