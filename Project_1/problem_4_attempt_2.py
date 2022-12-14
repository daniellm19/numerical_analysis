import numpy as np 
from numpy import sin, cos, sqrt, pi, linalg as LA
import matplotlib.pyplot as pst
from sympy.utilities.iterables import multiset_permutations

# All defined constants for the newton method
rho = 26570     # altitude of satellites [km]
c = 299792.458  # speed of light [km/s]
X = 0           # [km]
Y = 0           # [km]  # X, Y, and Z are the north pole (where the receiver is)
Z = 6370        # [km]
D = 0           # [s]   # This can be anything, just starts at 0

# Errors
err = 1e-8
toll_err = 1e-6
# An example location of satellites
corr_phi = [pi/8, pi/6, 3*pi/8, pi/4]       # Altitude
corr_theta = [-pi/4, pi/2, 2*pi/3, pi/6]    # Polar angle (azimuth)

def location(phi, theta):
    A = rho * sin(phi) * cos(theta) if abs(rho * sin(phi) * cos(theta)) > 1e-10 else 0
    B = rho * sin(phi) * sin(theta) if abs(rho * sin(phi) * sin(theta)) > 1e-10 else 0
    C = rho * cos(phi) if abs(rho * cos(phi)) > 1e-10 else 0
    distance = sqrt(pow((X - A), 2) + pow((Y - B), 2) + pow((Z - C), 2))
    t = distance / c

    ret_dict = {'A':A, 'B':B, 'C':C, 'distance':distance, 't':t}
    return ret_dict

def getABCt(incorr_phi):
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

def newtonmult(x0, tol, incorr_phi):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert r???? fyrir a?? F(x) og Jacobi fylki DF(x) s??u skilgreind annars sta??ar'''
    x=x0
    oldx=x+2*tol
    A, B, C, t = getABCt(incorr_phi)
    while LA.norm(x-oldx, np.inf)>tol:
        oldx=x
        s=-LA.solve(DF(x, A, B, C, t),F(x, A, B, C, t))
        x=x+s
    return(x)


def add_plus_min_err_to_list(err: int, a_list: list):
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
    all_perms = [[a + b for a, b in zip(all_err_perm[i], corr_phi)] for i in range(len(all_err_perm))]

    return all_perms

def main():
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    x0 = np.array([X,Y,Z,D]) #Initial guess for newtons method
    incorr_phis = add_plus_min_err_to_list(err, corr_phi)
    all_lenghts = []
    for incorr_phi in incorr_phis:
        x,y,z,d = newtonmult(x0, toll_err, incorr_phi)
        print('The error is:')
        print(sqrt(pow(x - X, 2) + pow(y - Y, 2) + pow(z - Z, 2)))
        all_lenghts.append(sqrt(pow(x - X, 2) + pow(y - Y, 2) + pow(z - Z, 2)))
    print(all_lenghts)
    
main()

if __name__ == "__main__":
    print(main)