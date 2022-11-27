import numpy as np 
from Positions import Positions
from numpy import sin, cos, sqrt, pi, linalg as LA
import matplotlib.pyplot as pst

if __name__ == "__main__":
    corr_phi = [pi/8, pi/6, 3*pi/8, pi/4]
    corr_theta = [-pi/4, pi/2, 2*pi/3, pi/6]
    incorr_phi = [pi/8 + 1e-8, pi/6 + 1e-8, 3*pi/8 - 1e-8, pi/4 - 1e-8]
    
    pos = Positions(np.array([0,0,6730,0]), 1e-8)
    pos.add_satellite_polar()
    print(pos.cartesian_calc(np.pi/3, np.pi/2))