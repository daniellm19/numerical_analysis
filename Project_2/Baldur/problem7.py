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
    m = 1
    y_0 = np.array([pi/3, 0, pi/6, 0])
    problem = input('Do you want to run problem 3 or 4? ')
    if problem == '4':
        y_0[0] = pi/2
        print('Running problem 4, enjoy!')
    else:
        print('Running problem 3, enjoy!')
    y, t, h = runge_kutta(y_0, n, T, L, m)
    print(y)
    position1, velocity1, position2, velocity2 = map(list, zip(*y))
    #ani_plot(t, position1, L, h)
    #plt.clf()
    plt.plot(t,y)
    plt.show()


main()
