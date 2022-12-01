import matplotlib.pyplot as plt
from math import sin, cos, pi
import numpy as np
from problem_2 import Problem_one_euler_method
import matplotlib.animation as animation

T = 20 # the interval of t [0, T]
n = 500 # the segments to break down the interval into
z_0 = [pi/12, 0] # inital values of pendulum in vector. Can be changed at convinience. [angular position, angular speed]
L = 2

bob_radius = 0.1
fig = plt.figure()
#fig.add_axes(3, 3, 3, 3)
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

plt.show()