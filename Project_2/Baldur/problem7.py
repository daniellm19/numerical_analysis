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
    
def animate_pendulum(angle1, velocity1, angle2, velocity2, L, T):
    bob_radius = 0.1
    
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    # Set the plot limits so that the pendulum has room to swing!
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)

    theta0 = pi/12

    def get_coords(th):
        """Return the (x, y) coordinates of the bob at angle th."""
        return L * np.sin(th), -L * np.cos(th)

    
    line, = ax.plot([0, -2*sin(angle1[0])], [0, 2*cos(velocity1[0])], lw=3, c='k') # initial pendulum position
    circle = ax.add_patch(plt.Circle(get_coords(angle1[0]), bob_radius,
                        fc='r', zorder=3))
    
    line_2, = ax.plot([0, -2*sin(angle2[0])], [0, 2*cos(velocity2[0])], lw=3, c='k') # initial pendulum 2 position
    circle_2 = ax.add_patch(plt.Circle(get_coords(angle2[0]), bob_radius,
                        fc='b', zorder=3))
    
    def animate(i):
        x_1 = 2 * sin(angle1[i])
        y_1 = 2 * -cos(angle1[i])
        x_2 = x_1 +  2 * sin(angle2[i])
        y_2 = y_1 + 2 * -cos(angle2[i])
        line.set_data([0, x_1], [0, y_1])
        circle.set_center((x_1, y_1))
        line_2.set_data([x_1, x_2], [y_1, y_2])
        circle_2.set_center((x_2, y_2))


    ani = animation.FuncAnimation(fig, animate, frames=1000, repeat=True, interval=T)
    plt.grid()
    plt.show()

def main():
    T = 20
    n = 500
    L = 2
    m = 1
    y_0 = np.array([pi+0.1, 0, pi, 0])
    y, t, h = runge_kutta(y_0, n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    animate_pendulum(angle1, velocity1, angle2, velocity2, L, T)


main()