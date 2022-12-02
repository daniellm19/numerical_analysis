import matplotlib.pyplot as plt
from numpy import sin, cos, pi
import numpy as np
import matplotlib.animation as animation

# Problem 2 starts here: 

def eulerstep(t, y, h):
    return  [y[0] + h*y_dot(t, y)[0], y[1] + h*y_dot(t, y)[1]]

def y_dot(t, y):
    g = 9.81 # force of gravity in m*s^-2
    L = 2 # lenght of pendulum in m
    return [y[1], -(g/L)*sin(y[0])]

def Problem_one_euler_method(T, n, y_0):
    """
    T: interval from [0, T]
    n: number of segmentations on interval [0, T]
    y_0: initial vector values for angular position and angular speed. [theta, theta_prime]
    
    This function uses the euler method to emulate the movement of pendumulum with differential equation from problem one (see equations in eulerstep and y_dot).
    Returns t, theta[0], theta[1] (t = time array, theta[0] = theta position array, theta[1] = theta speed array)
    """
    z = y_0
    h = T/n
    t = [i*h for i in range(0, n)]
    theta = [[], []] # array for the changing values of theta
    for i in range(0, n):
        z = eulerstep(t[i], z, h)
        theta[0].append(z[0])
        theta[1].append(z[1])

    return t, theta[0], theta[1]

# Problem 3 & 4 starts here

def animate_pendulum(x, y, h):
    
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-2.2, 2.2), ylim=(-2.2, 2.2))
    ax.grid()   
    
    line = ax.plot([], [], 'o-', c='blue', lw=1.5)
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

def plot(t,pos,vel):
    plt.figure(figsize=(8,4))
    plt.plot(t, pos, label="Pendulum's angle [rad]")
    plt.plot(t, vel, label = "Pendulum's angular velocity [rad/s]")
    plt.xlabel('Time [s]')
    plt.ylabel('Radians')
    plt.legend()
    plt.show()
    
def main():
    L = 2
    T = 20 # the interval of t [0, T]
    n = 500 # the segments to break down the interval into
    h = T/n
    
    problem = input("Do you want to run problem 3 or 4? (3/4): ")
    if problem == "3":
        y_0 = [pi/12, 0] # inital values of pendulum in vector for problem 3
    elif problem == "4":
        y_0 = [pi/2, 0] # inital values of pendulum in vector for problem 4
    else:
        print("wrong input, I´m shutting down!")
        quit()
        
    t, angle, velocity = Problem_one_euler_method(T, n, y_0) # get the estimation for eulers method on problem one differential equation
    x, y = L * sin(angle[:]), -L * cos(angle[:])
    
    animate_pendulum(x, y, h)
    plot(t, angle, velocity)
    
main()