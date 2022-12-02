import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos, pi


#Global constants
g = 9.81

def error(k):
    return pow(10,-k)

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
    T = 40
    n = 500
    L = 2
    m = 1
    k = [1,2,3,4,5]
    y0 = np.array([2*pi/3, 0, pi/6, 0])
    y, t, h = runge_kutta(y0, n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    diff_y = []
    counter = 1
    for num in k:
        y0 = np.array([2*pi/3+error(num), 0, pi/6 + error(num), 0])
        y_n, t_n, h_n = runge_kutta(y0, n, T, L, m)
        diff_y.append([y_n, t_n, h_n])
        angle1n, velocity1n, angle2n, velocity2n = map(list, zip(*y_n))
        plt.figure()
        plt.plot(t_n, angle1n)
        plt.plot(t, angle1)
        plt.grid()
        counter += 1  
        
    plt.show()
    #angle1, velocity1, angle2, velocity2 = map(list, zip(*y))

main()