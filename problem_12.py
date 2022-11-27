import numpy as np 
from numpy import sin, cos, sqrt, pi, mean, std, linalg as LA
import matplotlib.pyplot as plt
import random
from sympy.utilities.iterables import multiset_permutations
import json

#constants
rho = 26570 # kilometers
c = 299792.458 # speed of light [km/s]
X = 0
Y = 0
Z = 6370
err = 1e-8

def random_angles():
    '''Reads the random angles from random_angles.json'''
    with open('random_angles.json') as f:
        all_angles = json.load(f)
    rand_thetas, rand_phis = all_angles[0], all_angles[1]
    return rand_thetas, rand_phis

def location(phi: float, theta: float):
    A = rho * sin(phi) * cos(theta) if abs(rho * sin(phi) * cos(theta)) > 1e-10 else 0
    B = rho * sin(phi) * sin(theta) if abs(rho * sin(phi) * sin(theta)) > 1e-10 else 0
    C = rho * cos(phi) if abs(rho * cos(phi)) > 1e-10 else 0
    distance = sqrt(pow((X - A), 2) + pow((Y - B), 2) + pow((Z - C), 2))
    t = distance / c

    ret_dict = {'A':A, 'B':B, 'C':C, 'distance':distance, 't':t}
    return ret_dict


def getabct(corr_theta: list, corr_phi: list, incorr_angle: list, angle: str):
    A, B, C, t = [], [], [], []
    for i in range(len(corr_phi)):
        t.append(location(corr_phi[i], corr_theta[i])['t']) # Vector of time for each sat t[s] derived from correct values
        if angle == 'phi':
            values = location(incorr_angle[i], corr_theta[i])     # Derived from perceived values
            A.append(values['A'])                               # Vector of distances in plane A[km]
            B.append(values['B'])                               # Vector of distances in plane B[km]
            C.append(values['C'])                               # Vector of distances in plane C[km]
        else:
            values = location(corr_phi[i], incorr_angle[i])     # Derived from perceived values
            A.append(values['A'])                               # Vector of distances in plane A[km]
            B.append(values['B'])                               # Vector of distances in plane B[km]
            C.append(values['C'])                               # Vector of distances in plane C[km]
    return A, B, C, t

def F(x: list, A: list, B: list, C: list, t: list):
    funcs = []
    for i in range(4):
        funcs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(c,2) * pow((t[i]-x[3]), 2))
    return funcs

def DF(x: list, A: list, B: list, C: list, t: list):
    '''Creates each row of jacobi matrix independently(Hard coded) takes in a 4x1 vector
    of initial conditions and creates a matrix which is returned'''
    eq_list = [0,0,0,0]
    for i in range(0, len(eq_list)):
        eq_list[i] = [(2*(x[0]-A[i])), (2*(x[1]-B[i])), (2*(x[2]-C[i])), (2*pow(c,2)* (t[i]-x[3]))]
    return np.array(eq_list)

def newtonmult(x0: list , tol: int, theta: list, phi: list, incorr_angle: list, angle: str):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert ráð fyrir að F(x) og Jacobi fylki DF(x) séu skilgreind annars staðar'''
    x=x0
    oldx=x+2*tol
    A, B, C, t = getabct(theta, phi, incorr_angle, angle)
    while LA.norm(x-oldx, np.inf)>tol:
        #print(LA.norm(x-oldx, np.inf))
        oldx=x
        s=-LA.solve(DF(x, A, B, C, t), F(x, A, B, C, t))
        x=x+s
    return(x)

def distance_w_error(theta: list, phi: list, incorr_angle: list, angle: str):
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    x,y,z,d = newtonmult(x0, 1e-8, theta, phi, incorr_angle, angle)
    print(theta, phi, incorr_angle, angle)
    return (sqrt(pow(x - X, 2) + pow(y - Y, 2) + pow(z - Z, 2)))



def non_angle(a_list: list, wrong_value: float, how_many):
    b_list = []
    for i in range(0, how_many):
        b_list.append(wrong_value)

    for i in range(how_many, len(a_list)):
        b_list.append(a_list[i])
    return b_list

def main():
    replacement_theta = pi
    replacement_phi = 3*pi/16
    all_errors_wrong_thetas = []
    all_errors_wrong_phis = []
    rand_thetas, rand_phis = random_angles()

    for i in range(len(rand_thetas)):
        print(i)
        wrong_theta = non_angle(rand_thetas[i], replacement_theta, 3)
        max_error = distance_w_error(rand_thetas[i], rand_phis[i], wrong_theta, 'theta')
        all_errors_wrong_thetas.append(max_error)

    for i in range(0, 10):
        wrong_phi = non_angle(rand_phis[i], replacement_phi, 1)
        max_error = distance_w_error(rand_thetas[i], rand_phis[i], wrong_phi, 'phi')
        all_errors_wrong_phis.append(max_error)


    plt.clf()
    plt.subplot(121)   
    plt.boxplot(all_errors_wrong_thetas) 
    plt.ylabel('Perceived error [km]') 
    plt.subplot(122)
    plt.hist(all_errors_wrong_thetas)
    plt.xlabel('Perceived error [km]')
    plt.ylabel('No. of sattelite groups')
    plt.show()

    plt.clf()
    plt.subplot(121)   
    plt.boxplot(all_errors_wrong_phis) 
    plt.ylabel('Perceived error [km]') 
    plt.subplot(122)
    plt.hist(all_errors_wrong_phis)
    plt.xlabel('Perceived error [km]')
    plt.ylabel('No. of sattelite groups')
    plt.show()

import time
if __name__ == "__main__":
    main()



