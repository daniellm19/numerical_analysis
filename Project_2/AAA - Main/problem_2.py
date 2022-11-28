import matplotlib.pyplot as mp
from math import sin, pi
import numpy as np

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
    


