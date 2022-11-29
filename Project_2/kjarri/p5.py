import matplotlib.pyplot as plt
from numpy import sin, pi, cos, sqrt
import numpy as np


#constants:
g=9.81

def ydot(t:float, y: list, L: float):
    return np.array([y[1], -g/L*sin(y[0])])


def runge_kutta(x, n, T, L: float):
    h = T/n
    t = 0
    t_list = []
    y_list = []
    for _ in range(n):
        k1 = ydot(t, x, L)
        k2 = ydot(t + h/2, x + h/2 * k1, L)
        k3 = ydot(t + h/2, x + h/2 * k2, L)
        k4 = ydot(t + h, x + h * k3, L)
        y = np.array(x + h * (k1/6 + k2/3 + k3/3 + k4/6))
        y_list.append(y)
        x = y
        t += h
        t_list.append(t)


    return y_list, t_list, h

def point_coordenents(position: float, L: float):
    x = sin(position) * L
    y = -cos(position) * L
    return x, y

def ani_plot(t: list, positions: list, L: float, h: float):
    
    for position in positions:
        plt.clf()
        plt.grid()
        plt.xlim(-2.5,2.5)
        plt.ylim(-2.5,2.5)
        x, y = point_coordenents(position, L)
        plt.plot([0,x], [0,y])
        plt.scatter(x, y, s=100, c='r')
        plt.pause(h)
    plt.show()

def main():
    T = 20
    n = 500
    L = 2
    y_0 = np.array([pi/12,0])
    problem = input('Do you want to run problem 3 or 4? ')
    if problem == '4':
        y_0[0] = pi/2
        print('Running problem 4, enjoy!')
    else:
        print('Running problem 3, enjoy!')
    y, t, h = runge_kutta(y_0, n, T, L)
    position, velocity = map(list, zip(*y))
    ani_plot(t, position, L, h)
    plt.clf()
    plt.plot(t,y)
    plt.show()


main()
