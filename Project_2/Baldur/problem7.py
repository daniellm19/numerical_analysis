import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos, pi
<<<<<<< HEAD
=======
import matplotlib.animation as animation
>>>>>>> 8587107ac0ccd1be678e192378f535fe69695c6b

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
<<<<<<< HEAD

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
        plt.pause(h)
    plt.show()
=======
    
def animate_penduli(x_1, y_1, x_2, y_2, T):
    
    def animate(i):
        ax.clear()
        
        if i > len(x_1) - 1:
            print("Animation completed!\nThanks for playing")
            quit()
        # updating trajectory
        ax.plot([0, x_1[i]], [0, y_1[i]], lw='1', c='k')
        ax.plot(x_1[:i+1], y_1[:i+1], 'red', lw='0.5')
        
        # updating point location
        ax.scatter(x_1[i], y_1[i], c='red')
        
        # same but for second pendulum
        ax.plot([x_1[i], x_2[i]], [y_1[i], y_2[i]], lw='1', c='k')
        ax.plot(x_2[:i+1], y_2[:i+1], c='blue', lw='0.5')
        ax.scatter(x_2[i], y_2[i], c='blue')
            
        ax.set_xlim([-4, 4])
        ax.set_ylim([-4, 4])
        plt.grid()
        ax.set_xlabel('x')
        ax.set_xlabel('y')
        
    fig = plt.figure()
    ax = plt.axes()
   
    show_me_the_penduli = animation.FuncAnimation(fig, animate, frames=1000, repeat=False, interval=0.50)
    
    plt.grid()
    plt.show()
    
    
    
def get_coords(theta, L):
        """Return the (x, y) arrays of coordinates of the ball at angle theta."""
        x, y = [], []
        for i in theta:
            x.append(L * sin(i))
            y.append(L * -cos(i))
        return x, y

def second_pend_coords_offset(x_1, y_1, x_2, y_2):
    for i in range(len(x_1)):
        x_2[i] = x_1[i] + x_2[i]
        y_2[i] = y_1[i] + y_2[i]
    return x_2, y_2
>>>>>>> 8587107ac0ccd1be678e192378f535fe69695c6b

def main():
    T = 20
    n = 500
    L = 2
    m = 1
<<<<<<< HEAD
    y_0 = np.array([pi+0.1, 0, pi, 0])
    y, t, h = runge_kutta(y_0, n, T, L, m)
    print(y)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    ani_plot(t, angle1, angle2, L, h)
    plt.clf()



main()
=======
    y_0 = np.array([pi/3, 0, pi/6, 0])
    y, t, h = runge_kutta(y_0, n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    
    x_1, y_1 = get_coords(angle1, L)
    x_2, y_2 = get_coords(angle2, L, )
    x_2, y_2 = second_pend_coords_offset(x_1, y_1, x_2, y_2)
    
    animate_penduli(x_1, y_1, x_2, y_2, T)
    
    return

main()
>>>>>>> 8587107ac0ccd1be678e192378f535fe69695c6b
