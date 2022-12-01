import matplotlib.pyplot as plt
from math import sin, cos, pi
import numpy as np
import matplotlib.animation as animation

# Problem 2 starts here: 

def eulerstep(t, y, h):
    return  [y[0] + h*y_dot(t, y)[0], y[1] + h*y_dot(t, y)[1]]

def y_dot(t, y):
    g = 9.81 # force of gravity in m*s^-2
    L = 2 # lenght of pendulum in m
    return [y[1], -(g/L)*sin(y[0])]

def Problem_one_euler_method(T, n, z_0):
    """
    T: interval from [0, T]
    n: number of segmentations on interval [0, T]
    z_0: initial vector values for angular position and angular speed. [theta, theta_prime]
    
    This function uses the euler method to emulate the movement of pendumulum with differential equation from problem one (see equations in eulerstep and y_dot).
    Returns t, theta[0], theta[1] (t = time array, theta[0] = theta position array, theta[1] = theta speed array)
    """
    z = z_0
    h = T/n
    t = [i*h for i in range(0, n)]
    theta = [[], []] # array for the changing values of theta
    for i in range(0, n):
        z = eulerstep(t[i], z, h)
        theta[0].append(z[0])
        theta[1].append(z[1])

    return t, theta[0], theta[1]

# Problem 3 & 4 starts here

def animate_pendulum(z_0):
    
    L = 2
    T = 20 # the interval of t [0, T]
    n = 500 # the segments to break down the interval into
    bob_radius = 0.1
    
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    # Set the plot limits so that the pendulum has room to swing!
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)

    t, theta, theta_prime = Problem_one_euler_method(T, n, z_0) # get the estimation for eulers method on problem one differential equation
    theta0 = pi/12

    def get_coords(th):
        """Return the (x, y) coordinates of the bob at angle th."""
        return L * np.sin(th), -L * np.cos(th)

    line, = ax.plot([0, -2*sin(z_0[0])], [0, 2*cos(z_0[1])], lw=3, c='k') # initial pendulum position
    circle = ax.add_patch(plt.Circle(get_coords(theta0), bob_radius,
                        fc='r', zorder=3))
    def animate(i):
        x = 2 * sin(theta[i])
        y = 2 * -cos(theta[i])
        line.set_data([0, x], [0, y])
        circle.set_center((x, y))


    ani = animation.FuncAnimation(fig, animate, frames=1000, repeat=True, interval=T)
    plt.grid()
    plt.show()
    
def main():
    problem = input("Do you want to run problem 3 or 4? (3/4): ")
    if problem == "3":
        z_0 = [pi/12, 0] # inital values of pendulum in vector for problem 3
    elif problem == "4":
        z_0 = [pi/2, 0] # inital values of pendulum in vector for problem 4
    else:
        print("wrong input, I´m shutting down!")
        quit()
        
    animate_pendulum(z_0)
    
main()