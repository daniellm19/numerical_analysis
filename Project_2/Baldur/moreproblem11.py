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

def error_list(yn: list, y: list, n):
    ret_arr = []
    for i in range(n):
        ret_arr.append([abs(y[i][0]-yn[i][0]),abs(y[i][1]-yn[i][1]),abs(y[i][2]-yn[i][2]),abs(y[i][3]-yn[i][3])])
    return ret_arr

def print_the_error(error: list, t):
    fig = plt.figure(figsize=(12, 7))  
    plt.subplots_adjust(top=0.9, left=0.06, right=0.96, hspace=0.5, bottom=0.08)
    fig.suptitle("Difference in angle with starting error 10^(-k) with k ∈ {1,2,3,4,5}")
    counter = 1
    for element in error:
        plt.subplot(3, 2, counter)
        plt.grid()
        angle1, velocity1, angle2, velocity2 = map(list, zip(*element))
        plt.plot(t, angle1, label="Inner pendulum error")
        plt.plot(t, angle2, label="Outer pendulum error")
        plt.xlabel('Time [s]',fontsize=7)
        plt.ylabel('Radians', fontsize=7)
        plt.legend(loc=2, prop={'size': 6})
        plt.title(f"Starting with error 10^-({counter}) [rad]", fontsize=9)
        counter += 1 

def make_table(error: list):
    print("Error at T=40")
    print("Inner\t\t  Outer")
    for item in error:
        print(item[-1][0], item[-1][2])

def main():
    T = 40
    n = 500
    L = 2
    m = 1
    k = [1, 2, 3, 4, 5]
    y0 = np.array([2*pi/3, 0, pi/6, 0])
    y, t, h = runge_kutta(y0, n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    new_arr = []
    for num in k:
        y0 = np.array([2*pi/3+error(num), 0, pi/6 + error(num), 0])
        y_n, t_n, h_n = runge_kutta(y0, n, T, L, m)
        y_n_error = error_list(y_n, y, n)
        new_arr.append(y_n_error)
    print_the_error(new_arr, t)
    plt.show()  
    make_table(new_arr)
    #angle1, velocity1, angle2, velocity2 = map(list, zip(*y))

main()