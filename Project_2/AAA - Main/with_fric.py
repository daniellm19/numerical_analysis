import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, pi
import matplotlib.animation as animation
from collections import deque

def ydot(t:float, y: list, L: float, m:float, g:float, d:float):
    theta_1pp = ((m*L*pow(y[1],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m*g*sin(y[2])*cos(y[2] - y[0])) + m * L * pow(y[3],2) * sin(y[2] - y[0]) - (m+m) * g*sin(y[0]))/((m+m)*L - m*L*pow(cos(y[2] - y[0]),2)) - d * y[1]
    theta_2pp = (-1*(m*L*pow(y[3],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m+m)*(g*sin(y[0])*cos(y[2] - y[0])-L*pow(y[1], 2) * sin(y[2] - y[0])-g*sin(y[2])))/((m+m)*L - m*L*pow(cos(y[2] - y[0]), 2))
    return np.array([y[1], theta_1pp, y[3], theta_2pp])


def runge_kutta(x, n, T, L: float, m: float, g: float = 9.81, d: float = 0, E: float = 0):
    h = T/n
    t = 0
    t_list = []
    y_list = []
    initial_energy = energy(x, L, m, g)
    for _ in range(n):
        k1 = ydot(t, x, L, m, g, d)
        k2 = ydot(t + h/2, x + h/2 * k1, L, m, g, d)
        k3 = ydot(t + h/2, x + h/2 * k2, L, m, g, d)
        k4 = ydot(t + h, x + h * k3, L, m, g, d)
        y = np.array(x + h * (k1/6 + k2/3 + k3/3 + k4/6))
        y_list.append(y)
        x = y
        t += h
        t_list.append(t)
        if initial_energy - energy(x, L, m, g) >= E:
            t_list, y_th_list, x_th_list = get_throw(y_list, t_list, h, T, g)
            return y_list, t_list, h, y_th_list, x_th_list
            break
    return y_list, t_list, h

def energy(y, m, g, l):
    kin_e = (1/2) * (m + m) * pow(l, 2) * pow(y[1], 2) + (1/2) * m * (pow(l, 2) * pow(y[3], 2) + m * l * l * y[1] * y[3] * cos(y[0] - y[2]))
    pot_e = - (m + m) * g * l * cos(y[0]) - m * g * l * cos(y[3])
    
    return kin_e + pot_e

def get_throw(y_list, T_list, h, T, g):
    time = 0
    vy = - cos(y_list[-1][1])
    vx = sin(y_list[-1][1])
    x = vx * time
    y = vy * time - ((1/2) * g * pow(time, 2))
    x_th_list, y_th_list = [], []
    while abs(x) < 4 and abs(y) < 4:
        time += h
        x = vx * time
        y = vy * time - ((1/2) * g * pow(time, 2))
        T_list.append(time)
        y_th_list.append(y)
        x_th_list.append(x)
        print(x)
        print(y)
    return T_list, y_th_list, x_th_list
    


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
    return 

def plot(t,pos,vel, pendulum_name):
    plt.figure(figsize=(8,4))
    plt.plot(t, pos, label=f"{pendulum_name} angle [rad]")
    plt.plot(t, vel, label = f"{pendulum_name} angular velocity [rad/s]")
    plt.xlabel('Time [s]')
    plt.ylabel('Radians')
    plt.legend()
    plt.show()


def get_g():
    g = {'Sun': 274.13, 'Mercury': 3.7, 'Venus': 8.87, 'Earth': 9.81, 'Moon': 1.62, 'Mars': 3.71, 'Jupiter': 24.92, 'Saturn': 10.44, 'Uranus': 8.87, 'Neptune': 11.15, 'Pluto': 0.58}
    planet = input('Please enter a corresponding number for celestial body where the pendulum should be.\n' + 
    '0: Sun, 1: Mercury, 2: Venus, 3: Earth, 4: Moon, 5: Mars, 6: Jupiter, 7: Saturn, 8: Uranus, 9: Neptune\n' )
    try:
        print('You chose', list(g)[int(planet)], 'where the gravitational acceleration is', list(g.values())[int(planet)], 'm/s^2')
        return list(g.values())[int(planet)]
    except ValueError or IndexError:
        print('Please enter a number from 0-10')
        quit()

def get_E(m):
    metals = {'aluminum': [658, 910], 'lead': [327, 130], 'tin': [231, 210], 'silver': [960, 234], 'zinc': [419, 390]}
    metal = input("Please enter a corresponding number for metal that the pendulum's pivot point should be made of.\n" + 
    '0: Aluminum, 1: Lead, 2: Tin, 3: Silver, 4: Zinc\n' )
    try:
        print('You chose', list(metals)[int(metal)], 'with melting point of ', list(metals.values())[int(metal)][0], '°C and specific heat of', list(metals.values())[int(metal)][1], 'J/Kg°C')
    except ValueError or IndexError:
        print('Please enter a number from 0-11')
        quit()
    return list(metals.values())[int(metal)][1] * m * (list(metals.values())[int(metal)][0] - 25)


def main():
    T = 20
    n = 500
    L = 2 # Length in meters
    m = 1 # mass in kg
    m_pivot = 0.01 # Mass of the pivot point in Kg
    g = get_g() # Gravitational acceleration [m/s^2]
    d = 5 # Friction coefficient
    E = get_E(m_pivot) # How much energy is requiered to melt the pivot point
    y_0 = np.array([-pi+0.1, 0, -pi, 0])

    try:
        y, t, h, y_throw, x_throw = runge_kutta(y_0, n, T, L, m, g, d, E)
    except ValueError:
        y, t, h = runge_kutta(y_0, n, T, L, m, g, d, E)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    
    # get the x, y co-ordinates from angle positions
    x_1, y_1 = L * sin(angle1[:]), -L * cos(angle1[:])

    try:
        for i in range(len(x_throw)):
            np.append(x_1, x_throw[i])
            np.append(y_1, y_throw[i])

    except UnboundLocalError:
        pass

    x_2, y_2 = L * sin(angle2[:]) + x_1, -L * cos(angle2[:]) + y_1

    try:
        for i in range(len(x_throw)):
            np.append(x_2, L * sin(angle2[-1]) + x_1[-i])
            np.append(y_2, -L * cos(angle2[-1]) + y_1[-i])
    except UnboundLocalError:
        pass
    
    animate_penduli(x_1, y_1, x_2, y_2, n, h)




main()