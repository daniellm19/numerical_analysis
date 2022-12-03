import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera
matplotlib.use('Agg')
from numpy import sin,cos, pi
import matplotlib.animation as animation
import random

#Global constants
g = 9.81

def ydot(y: list, L: float, m:float):
    theta_1pp = ((m*L*pow(y[1],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m*g*sin(y[2])*cos(y[2] - y[0])) + m * L * pow(y[3],2) * sin(y[2] - y[0]) - (m+m) * g*sin(y[0]))/((m+m)*L - m*L*pow(cos(y[2] - y[0]),2))
    theta_2pp = (-1*(m*L*pow(y[3],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m+m)*(g*sin(y[0])*cos(y[2] - y[0])-L*pow(y[1], 2) * sin(y[2] - y[0])-g*sin(y[2])))/((m+m)*L - m*L*pow(cos(y[2] - y[0]), 2))
    return np.array([y[1], theta_1pp, y[3], theta_2pp])

def runge_kutta_plot(x, n, T, L: float, m, name, color, first_num, second_num):
    plt.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca()
    camera = Camera(fig)
    h = T/n
    t_list = [h*(i+1) for i in range(n)]
    y_list = []
    theta_1_list = []
    theta_2_list = []
    for i in range(n):
        k1 = ydot(x, L, m)
        k2 = ydot(x + h/2 * k1, L, m)
        k3 = ydot(x + h/2 * k2, L, m)
        k4 = ydot(x + h * k3, L, m)
        y = np.array(x + h * (k1/6 + k2/3 + k3/3 + k4/6))
        theta_1_list.append(y[0])
        theta_2_list.append(y[2])
        ax.plot(theta_1_list, theta_2_list, color=color)
        camera.snap()
        #plt.plot(y[0], y[2], markersize=3, marker='o')
        y_list.append(y)
        x = y

    plt.ylabel(r"$\theta_2$")
    plt.xlabel(r"$\theta_1$")
    theta_1_latex = r"$\theta_1(0)=$"+r"$\frac{\pi}{"+f"{first_num}"+r"}$"
    theta_2_latex = r"$\theta_2(0)=$"+r"$\frac{\pi}{"+f"{second_num}"+r"}$"
    plt.title("Parametrized curve in "+r"$\mathbb{R}^2$"+" for " + theta_1_latex + " and " + theta_2_latex)
    plt.legend([r"$\langle\theta_1(t), \theta_2(t)\rangle$"])
    anim = camera.animate(blit=False, interval=10)
    anim.save(f'{name}.mp4')


    return y_list, t_list

def main():
    i = 10
    t = []
    y = []
    T = 200
    n = 1500
    n2 = 10000
    L = 2
    m = 1

    y_11 = np.array([0, 0, 0, 7])
    runge_kutta_plot(y_11, n2, T, L, m, "y11", "green", 3, 3)


main()