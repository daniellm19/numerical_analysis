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
    errors_a1, errors_a2, errors_v1, errors_v2 = [], [], [], []
    for n in ns:
        y_bad_est, t, h = runge_kutta(y_0, n, T, L, m)
        a1, v1, a2, v2 = map(list, zip(*y_bad_est))
        fin_bad_a1, fin_bad_a2, fin_bad_v1, fin_bad_v2 = a1[-1], a2[-1], v1[-1], v2[-1]
        errors_a1.append(abs(fin_good_a1 -  fin_bad_a1))
        errors_a2.append(abs(fin_good_a2 -  fin_bad_a2))
        errors_v1.append(abs(fin_good_v1 -  fin_bad_v1))
        errors_v2.append(abs(fin_good_v2 -  fin_bad_v2))
    #np.norm(errors_a1, errors_a2, errors_v1, errors_v2)
    return np.linalg.norm(np.array([errors_a1, errors_a2, errors_v1, errors_v2]))
    #la1 = np.polyfit(errors_v1, ns, 4)
    #print(la1)
    #la2 = np.linalg.lstsq(np.vstack([ns, np.ones(len(ns))]).T, errors_a2, rcond=None)[0]
    #lv1 = np.linalg.lstsq(np.vstack([ns, np.ones(len(ns))]).T, errors_v1, rcond=None)[0]
    #lv2 = np.linalg.lstsq(np.vstack([ns, np.ones(len(ns))]).T, errors_v2, rcond=None)[0]
    #return la1, la2, lv1, lv2



def main():
    T = 20
    good_n = 35000
    ns = np.array([100,200,400,800,1600,3200,6400])
    L = 2
    m = 1
    all_ini_vals = random_theta()[0:5]
    all_errors_theta_1, all_errors_theta_2, all_errors_thetap_1, all_errors_thetap_2 = [], [], [], []
    i=0
    for ini_val in all_ini_vals:
        print(i)
        ea1 = get_error(ini_val, good_n, ns, T, L, m)
        all_errors_theta_1.append(ea1)
        #all_errors_theta_2.append(ea2)
        #all_errors_thetap_1.append(ev1)
        #all_errors_thetap_2.append(ev2)
        i += 1
    print(all_errors_theta_1)
    plt.subplot(221)
    for ea1o in all_errors_theta_1:
        plt.plot(ns, ea1o)
    plt.ylabel('Error in angle θ1 [rad]')
    plt.xlabel('n')
    #plt.xscale('log', base=2)
    #plt.yscale('log', base=4)
    plt.subplot(222)
    for ea2 in all_errors_theta_2:
        plt.plot(ea2)
    plt.ylabel('Error in angle θ2 [rad]')
    plt.xlabel('n')
    plt.xscale('log', base=4)
    plt.yscale('log', base=4)
    plt.subplot(223)
    for ev1 in all_errors_thetap_1:
        plt.plot(ev1)
    plt.ylabel("Error in angular velocity θ'1 [rad/s]")
    plt.xlabel('n')
    plt.xscale('log', base=4)
    plt.yscale('log', base=4)
    plt.subplot(224)
    for ev2 in all_errors_thetap_2:
        plt.plot(ev2)
    plt.ylabel("Error in angular velocity θ'2 [rad/s]")
    plt.xlabel('n')
    plt.xscale('log', base=4)
    plt.yscale('log', base=4)

    plt.show()



main()