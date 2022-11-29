import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos, pi

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

def point_coordenents(angle1: float, angle2: float, L: float):
    x1 = sin(angle1) * L
    y1 = -cos(angle1) * L
    x2 = L * sin(angle2) + x1
    y2 = -cos(angle2) * L + y1
    return x1, x2, y1, y2

def ani_plot(t: list, angle1: list, angle2: list, L: float, h: float):
    plt.grid()
    for i in range(len(angle1)):
        plt.clf()
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        x1, x2, y1, y2 = point_coordenents(angle1[i], angle2[i], L)
        plt.plot([0,x1], [0,y1])
        plt.scatter(x1, y1, s=10, c='r')
        plt.scatter(x2, y2, s=10, c='b')
        plt.plot([x1,x2], [y1,y2])
        plt.pause(h/26)
    plt.show()

def main():
    T = 20
    n = 500
    L = 2
    m = 1
    y_0 = np.array([pi+0.1, 0, pi, 0])
    y, t, h = runge_kutta(y_0, n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    ani_plot(t, angle1, angle2, L, h)
    plt.clf()



main()
