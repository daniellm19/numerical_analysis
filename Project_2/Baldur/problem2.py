import numpy as np
from numpy import sin


#Global constants
g = 9.81 #Gravitational constant

def eulerstep(t,x,h):
    return y=x+h*y_prime(t,x)

def y_prime(t,y):
    L = 2
    return [y[1],-g/L*sin(y[0]]

def runner(T, init_theta, init_theta_prime, n):
    h = T/n
    y = [init_theta]
    t = [0]
    for i in range(n):

def func(theta):
    L = 2
    return g/l*sin(theta) 
    



