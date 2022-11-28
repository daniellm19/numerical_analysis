import matplotlib.pyplot as mp
from math import sin, pi
import numpy as np
from problem_2 import eulerstep

T = 20 # the interval of t [0, T]
n = 5000 # the segments to break down the interval into

def main():
    """
    This program defines initial theta and theta_prime z_0 (inital angular position and angular speed) and plots the theta 
    and theta_prime over time t.
    """
    z_0 = [pi/2, 0] # inital values of pendulum in vector. Can be changed at convinience. [angular position, angular speed]
    z = z_0
    h = T/n
    t = [i*h for i in range(0, n)]
    for i in range(0, n):
        z = eulerstep(t[i], z, h)
        theta[0].append(z[0])
        theta[1].append(z[1])
        
    mp.plot(t, theta[0])
    mp.plot(t, theta[1]) 
    mp.show()   
    
    
main()