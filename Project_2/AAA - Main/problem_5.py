import matplotlib.pyplot as plt
from numpy import sin, pi, cos, sqrt
import numpy as np
import matplotlib.animation as animation

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

def animate_pendulum(x, y, h):
    
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-2.2, 2.2), ylim=(-2.2, 2.2))
    ax.grid()   
    
    line, = ax.plot([], [], 'o-', c='blue', lw=1.5)
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
  
    def animate(i):
        xLine = [0, x[i]]
        yLine = [0, y[i]]
        
        line.set_data(xLine, yLine)
        time_text.set_text(f"time = {i*h:.1f}s")
        return line, time_text
   
    ani = animation.FuncAnimation(
        fig, animate, len(x), interval=h*1000, blit=True, repeat=False)
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
    angle, velocity = map(list, zip(*y))
    x, y = L * sin(angle[:]), -L * cos(angle[:])
    animate_pendulum(x, y, h)
    
    plt.clf()
    plt.plot(t,y)
    plt.show()


main()