import numpy as np 
from numpy import sin, cos, sqrt, pi, mean, std, linalg as LA
import matplotlib.pyplot as plt
import random
from sympy.utilities.iterables import multiset_permutations

#constants
rho = 26570 # kilometers
c = 299792.458 # speed of light [km/s]
X = 0
Y = 0
Z = 6370

def random_angles():
    '''Generates 200 random positions for four sattelites each'''
    rand_phis = []
    rand_thetas = []
    for i in range(0,100):
        rand_phi = []
        rand_theta = []
        for j in range(0,4):
            rand_phi.append(random.uniform(0.0, pi/2))
            rand_theta.append(random.uniform(0.0, 2*pi))
        rand_thetas.append(rand_theta)
        rand_phis.append(rand_phi)
    return rand_thetas, rand_phis

def location(phi, theta):
    A = rho * sin(phi) * cos(theta) if abs(rho * sin(phi) * cos(theta)) > 1e-10 else 0
    B = rho * sin(phi) * sin(theta) if abs(rho * sin(phi) * sin(theta)) > 1e-10 else 0
    C = rho * cos(phi) if abs(rho * cos(phi)) > 1e-10 else 0
    distance = sqrt(pow((X - A), 2) + pow((Y - B), 2) + pow((Z - C), 2))
    t = distance / c

    ret_dict = {'A':A, 'B':B, 'C':C, 'distance':distance, 't':t}
    return ret_dict



err = 1e-8

def get_incorr_phis(err: int, a_list: list):
    """Takes in a error and the list of values to add the combinations of errors to"""
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


'''
incorr_phis = [[pi/8 + err, pi/6 + err, 3*pi/8 + err, pi/4 + err], [pi/8 + err, pi/6 + err, 3*pi/8 + err, pi/4 - err],
             [pi/8 + err, pi/6 + err, 3*pi/8 - err, pi/4 + err], [pi/8 + err, pi/6 - err, 3*pi/8 + err, pi/4 + err], 
             [pi/8 - err, pi/6 + err, 3*pi/8 + err, pi/4 + err], [pi/8 - err, pi/6 - err, 3*pi/8 - err, pi/4 - err],
             [pi/8 - err, pi/6 - err, 3*pi/8 + err, pi/4 + err], [pi/8 + err, pi/6 - err, 3*pi/8 - err, pi/4 + err],
             [pi/8 + err, pi/6 + err, 3*pi/8 - err, pi/4 - err], [pi/8 - err, pi/6 + err, 3*pi/8 + err, pi/4 - err],
             [pi/8 + err, pi/6 - err, 3*pi/8 + err, pi/4 - err], [pi/8 - err, pi/6 + err, 3*pi/8 - err, pi/4 + err],
             [pi/8 + err, pi/6 - err, 3*pi/8 - err, pi/4 - err], [pi/8 - err, pi/6 + err, 3*pi/8 - err, pi/4 - err],
             [pi/8 - err, pi/6 - err, 3*pi/8 + err, pi/4 - err], [pi/8 - err, pi/6 - err, 3*pi/8 - err, pi/4 + err]]
'''

def getABCt(corr_theta, corr_phi, incorr_phi):
    A, B, C, t = [], [], [], []
    for i in range(len(corr_phi)):
        t.append(location(corr_phi[i], corr_theta[i])['t']) #Vector of time for each sat t[s] derived from correct values
        values = location(incorr_phi[i], corr_theta[i]) #Derived from prerceived values
        A.append(values['A']) #Vector of distances in plane A[km]
        B.append(values['B']) #Vector of distances in plane B[km]
        C.append(values['C']) #Vector of distances in plane C[km]
    return A, B, C, t

def F(x, A, B, C, t):
    funcs = []
    for i in range(4):
        funcs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(c,2) * pow((t[i]-x[3]), 2))
    return funcs

def DF(x, A, B, C, t):
    '''Creates each row of jacobi matrix independently(Hard coded) takes in a 4x1 vector
    of initial conditions and creates a matrix which is returned'''
    l1 = [(2*(x[0]-A[0])), (2*(x[1]-B[0])), (2*(x[2]-C[0])), (2*pow(c,2)* (t[0]-x[3]))]
    l2 = [(2*(x[0]-A[1])), (2*(x[1]-B[1])), (2*(x[2]-C[1])), (2*pow(c,2)* (t[1]-x[3]))]
    l3 = [(2*(x[0]-A[2])), (2*(x[1]-B[2])), (2*(x[2]-C[2])), (2*pow(c,2)* (t[2]-x[3]))]
    l4 = [(2*(x[0]-A[3])), (2*(x[1]-B[3])), (2*(x[2]-C[3])), (2*pow(c,2)* (t[3]-x[3]))]
    return np.array([l1,l2,l3,l4])

def newtonmult(x0, tol, theta, phi, incorr_phi):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert ráð fyrir að F(x) og Jacobi fylki DF(x) séu skilgreind annars staðar'''
    x=x0
    oldx=x+2*tol
    A, B, C, t = getABCt(theta, phi, incorr_phi,)
    while LA.norm(x-oldx, np.inf)>tol:
        oldx=x
        s=-LA.solve(DF(x, A, B, C, t),F(x, A, B, C, t))
        x=x+s
    return(x)

def distance_w_error(theta, phi):
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    incorr_phis = get_incorr_phis(err, phi)
    all_lenghts = []
    for incorr_phi in incorr_phis:
        x,y,z,d = newtonmult(x0, 1e-8, theta, phi, incorr_phi)
        all_lenghts.append(sqrt(pow(x - X, 2) + pow(y - Y, 2) + pow(z - Z, 2)))
    return all_lenghts

def main():
    all_errors = []
    rand_thetas, rand_phis = random_angles()
    for i in range(len(rand_phis)):
        max_error = max(distance_w_error(rand_thetas[i], rand_phis[i]))
        all_errors.append(max_error)
    print('Max error in distance:', max(all_errors))
    print('Min error in distance:', min(all_errors))
    print('Average:', mean(all_errors))
    print('Standard deviation:', std(all_errors))
    fig = plt.figure(figsize =(10, 7))    
    plt.boxplot(all_errors)  
    plt.show()

import time
if __name__ == "__main__":
    main()



