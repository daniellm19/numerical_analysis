import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos, pi
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
    t_list = []
    y_list = []
    for _ in range(n):
        k1 = ydot(t, x, L, m)
        k2 = ydot(t + h/2, x + h/2 * k1, L, m)
        k3 = ydot(t + h/2, x + h/2 * k2, L, m)
        k4 = ydot(t + h, x + h * k3, L, m)
        y = np.array(x + h * (k1/6 + k2/3 + k3/3 + k4/6))
        y_list.append(y)
        x = y
        t += h
        t_list.append(t)


    return y_list, t_list, h

def random_theta():
    '''Reads the random angles from random_ini_values.json'''
    with open('random_ini_values.json') as f:
        all_ini_vals = json.load(f)
    return all_ini_vals
    
def get_error(ini_val, good_n, ns, T, L, m):
    y_0 = np.array(ini_val)
    y_good_est, t, h = runge_kutta(y_0, good_n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y_good_est))
    fin_good_a1, fin_good_a2, fin_good_v1, fin_good_v2 = angle1[-1], angle2[-1], velocity1[-1], velocity2[-1]
    errors_norm = []
    for n in ns:
        y_bad_est, t, h = runge_kutta(y_0, n, T, L, m)
        a1, v1, a2, v2 = map(list, zip(*y_bad_est))
        #print(np.array([abs(fin_good_a1 - a1[-1]), abs(fin_good_a2 - a2[-1]), abs(fin_good_v1 - v1[-1]), abs(fin_good_v2 - v2[-1])]))
        errors_norm.append(np.linalg.norm(np.array([abs(fin_good_a1 - a1[-1]), abs(fin_good_a2 - a2[-1]), abs(fin_good_v1 - v1[-1]), abs(fin_good_v2 - v2[-1])])))
    print(errors_norm)
    return errors_norm




def main():
    T = 20
    good_n = 35000
    ns = np.array([100,200,400,800,1600,3200,6400])
    L = 2
    m = 1
    all_ini_vals = random_theta()[0:2]
    all_errors_norm = []
    i=0
    for ini_val in all_ini_vals:
        print(i)
        all_errors_norm.append(get_error(ini_val, good_n, ns, T, L, m))

        i += 1

    for e in all_errors_norm:
        plt.plot(ns, e)
    plt.ylabel('Error in angle θ1 [rad]')
    plt.xlabel('n')
    plt.xscale('log', base=2)
    plt.yscale('log', base=4)
    plt.show()



main()