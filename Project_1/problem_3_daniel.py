import numpy as np 
from Positions import Positions
from numpy import sin, cos, sqrt, pi, linalg as LA
import matplotlib.pyplot as pst

if __name__ == "__main__":
    corr_phi = [np.pi/8, np.pi/6, 3*np.pi/8, np.pi/4]
    corr_theta = [-np.pi/4, np.pi/2, 2*np.pi/3, np.pi/6]
    incorr_phi = [np.pi/8 + 1e-8, np.pi/6 + 1e-8, 3*np.pi/8 - 1e-8, np.pi/4 - 1e-8]
    
    pos = Positions(np.array([0,0,6730,0]), toll=1e-8)
    A, B, C, t = [], [], [], []
    for i in range(len(corr_phi)):
        t = pos.cartesian_calc(corr_phi[i], corr_theta[i])['t'] #Vector of time for each sat t[s] derived from correct values
        abc_values = pos.cartesian_calc(incorr_phi[i], corr_theta[i]) #Derived from perrceived values
        pos.add_satellite_cart(abc_values['A'], abc_values['B'], abc_values['C'], t)

    x,y,z,d = pos.find_intersection()
    print(f"Values: x = {x} ,  y = {y} , z = {z} , d =  {d}")
    print('The error is: ', end="")
    print(np.sqrt(pow(x - pos.x0[0], 2) + pow(y - pos.x0[1], 2) + pow(z - pos.x0[2], 2)), 'km')