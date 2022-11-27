import numpy as np 
from Positions import Positions
from numpy import sin, cos, sqrt, pi, linalg as LA
import matplotlib.pyplot as pst

def main():
    # Problem 1
    print("Problem 1 solution: \n\n")
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    pos = Positions(pos_guess=x0, toll=1e-8)
    A = [15600, 18760, 17610, 19170]         # Vector of distances in plane A[km]
    B = [7540, 2750, 14630, 610]             # Vector of distances in plane B[km]
    C = [20140, 18610, 13480, 18390]         # Vector of distances in plane C[km]
    t = [0.07074, 0.07220, 0.07690, 0.07242] # Vector of time for each sat t[s]
    for i in range(0,4):
        pos.add_satellite_cart(A[i], B[i], C[i], t[i])
    x,y,z,d = pos.find_intersection()
    print("Final answers: \n")
    print("x is {:.4f} \t".format(x))
    print("y is {:.4f} \t".format(y))
    print("z is {:.2f} \t".format(z))
    print("d is {:.8f} \n".format(d))

    # Problem 2
    print("Problem 2 solution: \n\n")
    pos = Positions(np.array([0,0,6730,0]), 1e-8)
    print(pos.cartesian_calc(np.pi/3, np.pi/2))

    # Problem 3
    print("Problem 3 solution: \n\n")
    corr_phi = [pi/8, pi/6, 3*pi/8, pi/4]
    corr_theta = [-pi/4, pi/2, 2*pi/3, pi/6]
    incorr_phi = [pi/8 + 1e-8, pi/6 + 1e-8, 3*pi/8 - 1e-8, pi/4 - 1e-8]
    
    pos = Positions(np.array([0,0,6730,0]), 1e-8)
    pos.calculate_multiple_cartesian(corr_theta, corr_phi)
    print(pos.cartesian_calc(np.pi/3, np.pi/2))

    return  
main()