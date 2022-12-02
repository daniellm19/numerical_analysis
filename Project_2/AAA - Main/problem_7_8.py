import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
import matplotlib.animation as animation
from collections import deque

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

def animate_penduli(x_1, y_1, x_2, y_2, n, h):
    
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-4.2, 4.2), ylim=(-4.2, 4.2))
    ax.grid()   
    
    line_1, = ax.plot([], [], 'o-', c='blue', lw=1.5)
    line_2, = ax.plot([], [], 'o-', c='red', lw=1.5)
    trace, = ax.plot([], [], '.-', c='red', lw=0.5, ms=1)
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    trajectory_x, trajectory_y = deque(maxlen=n), deque(maxlen=n)
    
    def animate(i):
        x1 = [0, x_1[i]]
        y1 = [0, y_1[i]]

        x2 = [x_1[i], x_2[i]]
        y2 = [y_1[i], y_2[i]]
        
        if i == 0:
            trajectory_x.clear()
            trajectory_y.clear()

        trajectory_x.appendleft(x2[1])
        trajectory_y.appendleft(y2[1])

        
        line_1.set_data(x1, y1)
        line_2.set_data(x2, y2)
        trace.set_data(trajectory_x, trajectory_y)
        time_text.set_text(f"time = {i*h:.1f}s")
        return line_1, line_2, trace, time_text
   
    ani = animation.FuncAnimation(
        fig, animate, len(x_1), interval=h*1000, blit=True, repeat=False)
    plt.show()

def plot(t,pos,vel, pendulum_name):
    plt.figure(figsize=(8,4))
    plt.plot(t, pos, label=f"{pendulum_name} angle [rad]")
    plt.plot(t, vel, label = f"{pendulum_name} angular velocity [rad/s]")
    plt.xlabel('Time [s]')
    plt.ylabel('Radians')
    plt.legend()
    plt.show()

def main():
    T = 20
    n = 500
    L = 2
    m = 1
    y_0 = np.array([pi/3, 0, pi/6, 0])
    inital_values = input("\n1 ->   θ_1 = π/3 , θ'_1 = 0 , θ_1 = π/6 , θ'_1 = 0 . (These are the inital values for problem 7)\n"
        + "\nHere are some other random initial values:"
        + "\n2 ->   θ_1 = π , θ'_1 = 0 , θ_2 = π , θ'_2 = 4 "
        + "\n3 ->   θ_1 = π/12 , θ'_1 = 3 , θ_2 = π/6 , θ'_2 = -3 "
        + "\n4 ->   θ_1 = 0 , θ'_1 = 0 , θ_2 = 2π/3 , θ'_2 = 8 "
        
        "\n\nChoose inital values for both penduli: (1/2/3/4) " )
    if inital_values == '1':
        pass
    elif inital_values == '2':
        y_0 = np.array([pi, 0, pi, 4])
    elif inital_values == '3':
        y_0 = np.array([pi/12, 3, pi/6, -3])
    elif inital_values == '4':
        y_0 = np.array([0, 0, 2*pi/3, 8]) 
        
    y, t, h = runge_kutta(y_0, n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    
    # get the x, y co-ordinates from angle positions
    x_1, y_1 = L * sin(angle1[:]), -L * cos(angle1[:])
    x_2, y_2 = L * sin(angle2[:]) + x_1, -L * cos(angle2[:]) + y_1
    
    animate_penduli(x_1, y_1, x_2, y_2, n, h)
    plot(t, angle1, velocity1, "Pendulum 1 (blue pendulum)")
    plot(t, angle2, velocity2, "Pendulum 2 (red pendulum)")
    return

main()