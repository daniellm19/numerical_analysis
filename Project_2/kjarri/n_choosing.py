import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos, pi
import matplotlib.animation as animation

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
    

def main():
    T = 20
    #n = 500
    L = 2
    m = 1
    angle_1_final = []
    angle_2_final = []
    nn = []
    for n in range(10,50):
        y_0 = np.array([pi/3, 0, pi/6, 0])
        y, t, h = runge_kutta(y_0, n*1000, T, L, m)
        angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
        angle_1_final.append(angle1[-1])
        angle_2_final.append(angle2[-1])
        nn.append(n*1000)
        print(n)


    plt.subplot(121)
    plt.plot(nn, angle_1_final)
    plt.ylabel('Estimated angle θ1 [rad]')
    plt.xlabel('Value for n')
    plt.subplot(122)
    plt.plot(nn, angle_2_final)
    plt.ylabel('Estimated angle θ2 [rad]')
    plt.xlabel('Value for n')
    plt.show()




main()