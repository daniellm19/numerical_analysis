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

def random_angles(pos: int, sat_amount: int):
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

def get_incorr_phis(err: int, a_list: list):
    """Takes in a error and the list of values to add the combinations of errors to.
    Used for testing how positional errors affect the accuracy of the satellites"""
    all_err_perm = []
    err_list = [err] * len(a_list)
    for i in range(0,len(err_list)+1):
        err_copy = err_list.copy()
        for j in range(0,i):
            err_copy[j] = -err_list[j]
        multiset_perm = multiset_permutations(err_copy)
        for i in multiset_perm:
            all_err_perm.append(i)

    # Adds all of the error permutations to the a_list
    all_perms = [[a + b for a, b in zip(all_err_perm[i], a_list)] for i in range(len(all_err_perm))]

    return all_perms

def getabct(corr_theta: list, corr_phi: list, incorr_phi: list):
    A, B, C, t = [], [], [], []
    for i in range(len(corr_phi)):
        t.append(location(corr_phi[i], corr_theta[i])['t']) # Vector of time for each sat t[s] derived from correct values
        values = location(incorr_phi[i], corr_theta[i])     # Derived from perceived values
        A.append(values['A'])                               # Vector of distances in plane A[km]
        B.append(values['B'])                               # Vector of distances in plane B[km]
        C.append(values['C'])                               # Vector of distances in plane C[km]
    return A, B, C, t

def F(x: list, A: list, B: list, C: list, t: list):
    eqs = []
    for i in range(4):
        eqs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(c,2) * pow((t[i]-x[3]), 2))
    return eqs

def DF(x: list, A: list, B: list, C: list, t: list):
    '''Creates each row of Jacobi matrix independently (Hard coded) takes in a 4x1 vector
    of initial conditions and creates a matrix which is returned'''
    eq_list = [0,0,0,0]
    for i in range(0, len(eq_list)):
        eq_list[i] = [(2*(x[0]-A[i])), (2*(x[1]-B[i])), (2*(x[2]-C[i])), (2*pow(c,2)* (t[i]-x[3]))]
    return np.array(eq_list)

def newtonmult(x0: list , tol: int, theta: list, phi: list, incorr_phi: list):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert r???? fyrir a?? F(x) og Jacobi fylki DF(x) s??u skilgreind annars sta??ar'''
    x=x0
    oldx=x+2*tol
    A, B, C, t = getabct(theta, phi, incorr_phi)
    while LA.norm(x-oldx, np.inf)>tol:
        oldx=x
        s=-LA.solve(DF(x, A, B, C, t), F(x, A, B, C, t))
        x=x+s
    return(x)

def distance_w_error(theta: list, phi: list):
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    incorr_phis = get_incorr_phis(err, phi)
    all_lenghts = []
    for incorr_phi in incorr_phis:
        x,y,z,d = newtonmult(x0, 1e-8, theta, phi, incorr_phi)
        all_lenghts.append(sqrt(pow(x - X, 2) + pow(y - Y, 2) + pow(z - Z, 2)))
    return all_lenghts

def errors_no_outliers(a_list):
    pass
    no = []
    for err in a_list:
        if err < 0.00225:
            no.append(err)
    print(len(no))
    return no

def main():
    all_errors = []
    rand_thetas, rand_phis = random_angles(100,4)
    for i in range(len(rand_phis)):
        max_error = max(distance_w_error(rand_thetas[i], rand_phis[i]))
        all_errors.append(max_error)
    print(all_errors)
    print('Max error in distance:', max(all_errors))
    print('Min error in distance:', min(all_errors))
    print('Average:', mean(all_errors))
    print('Standard deviation:', std(all_errors))
    no_outliers = errors_no_outliers(all_errors)
    plt.clf()
    fig = plt.figure(figsize =(10, 7))    
    plt.boxplot(all_errors) 
    plt.ylabel('Perceived error [km]') 
    plt.show()
    plt.hist(all_errors)
    plt.xlabel('Perceived error [km]')
    plt.ylabel('No. of sattelite groups')
    plt.show()
    plt.hist(no_outliers)
    plt.xlabel('Perceived error [km]')
    plt.ylabel('No. of sattelite groups')
    plt.show()

import time
if __name__ == "__main__":
    main()



