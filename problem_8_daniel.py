import numpy as np
from numpy import sin, cos, sqrt, pi, mean, std, linalg as LA
import matplotlib.pyplot as plt
import random
from sympy.utilities.iterables import multiset_permutations


# Constants
rho = 26570 # kilometers
c = 299792.458 # speed of light [km/s]
X = 0
Y = 0
Z = 6370
err = 1e-8

def random_angles(pos: int, rows: int):
    '''Generates {pos} random positions for {rows} amount of satellites each'''
    rand_phis = []
    rand_thetas = []
    for i in range(0, pos):
        rand_phi = []
        rand_theta = []
        for j in range(0, rows):
            rand_phi.append(random.uniform(0.0, pi/2))
            rand_theta.append(random.uniform(0.0, 2*pi))
        rand_thetas.append(rand_theta)
        rand_phis.append(rand_phi)
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

def getABCt(corr_theta: list, corr_phi: list, incorr_phi: list):
    A, B, C, t = [], [], [], []
    for i in range(len(corr_phi)):
        t.append(location(corr_phi[i], corr_theta[i])['t']) #Vector of time for each sat t[s] derived from correct values
        values = location(incorr_phi[i], corr_theta[i]) #Derived from prerceived values
        A.append(values['A']) #Vector of distances in plane A[km]
        B.append(values['B']) #Vector of distances in plane B[km]
        C.append(values['C']) #Vector of distances in plane C[km]
    return A, B, C, t

def F(x: list, A: list, B: list, C: list, t: list, rows: int):
    '''Takes in a 5x1 vector as input for initial value. Creates a list and appends each equation 
    with relevant variable data and returns'''
    funcs = []
    for i in range(rows):
        funcs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(c,2) * pow((t[i]-x[3]), 2))
    return funcs

def DF(x: list, A: list, B: list, C: list, t: list, rows: int):
    '''Creates each row of jacobi matrix independently(Hard coded) takes in a 5x1 vector
    of initial conditions and creates a matrix which is returned'''
    eq_list = [0]*rows
    for i in range(0, len(eq_list)):
        eq_list[i] = [(2*(x[0]-A[i])), (2*(x[1]-B[i])), (2*(x[2]-C[i])), (2*pow(c,2)* (t[i]-x[3]))]
    return np.array(eq_list)

def newtonmult(x0,tol):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert ráð fyrir að F(x) og Jacobi fylki DF(x) séu skilgreind annars staðar'''
    x=x0
    oldx=x+2*tol
    while LA.norm(x-oldx, np.inf)>tol:
        oldx=x
        s=-LA.solve(DF(x),F(x))
        x=x+s
    return(x)

def newtonmult(x0: list , tol: int, theta: list, phi: list, incorr_phi: list):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert ráð fyrir að F(x) og Jacobi fylki DF(x) séu skilgreind annars staðar'''
    x=x0
    oldx=x+2*tol
    A, B, C, t = getABCt(theta, phi, incorr_phi)
    while LA.norm(x-oldx, np.inf)>tol:
        oldx=x
        s=-LA.solve(DF(x, A, B, C, t), F(x, A, B, C, t))
        x=x+s
    return(x)

def newton_gauss(x0: list, tol: int, theta: list, phi: list, incorr_phi: list, rows: int):
    """Implements the Newton-Gauss method, n are the amount of equations (rows) and m the amount of variables (columns)"""
    x = x0
    oldx=x+2*tol
    A, B, C, t = getABCt(theta, phi, incorr_phi)
    while LA.norm(x-oldx, np.inf)>tol:
        oldx = x
        left_side = np.dot(np.transpose(DF(x, A, B, C, t, rows)), DF(x, A, B, C, t, rows))
        right_side = np.dot((np.transpose(DF(x, A, B, C, t, rows))), F(x, A, B, C, t, rows))
        s = -LA.solve(left_side, right_side)
        x = x+s
    return(x)

def distance_w_error(theta: list, phi: list, err: int, sat_amount: int):
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    incorr_phis = get_incorr_phis(err, phi)
    all_lenghts = []
    for incorr_phi in incorr_phis:
        x,y,z,d = newton_gauss(x0, err, theta, phi, incorr_phi, sat_amount)
        all_lenghts.append(sqrt(pow(x - X, 2) + pow(y - Y, 2) + pow(z - Z, 2)))
    return all_lenghts

def main():
    sat_pos_amount = 100 
    sat_amount = 7
    all_all_errors = []
    max_errors = []
    min_errors = []
    mean_errors = []
    for sat_amount in range(5, sat_amount+1):
        all_errors = []

        rand_thetas, rand_phis = random_angles(sat_pos_amount, sat_amount)
        for i in range(len(rand_phis)):
            max_error = max(distance_w_error(rand_thetas[i], rand_phis[i], err, sat_amount))
            all_errors.append(max_error)

        max_error = max(all_errors)
        min_error = min(all_errors)
        mean_error = mean(all_errors)
        print('Max error in distance:', max_error)
        print('Min error in distance:', min_error)
        print('Mean error in distance:', mean_error)
        print(f'Standard deviation: {std(all_errors)}\n')

        fig, ax = plt.subplots()
        ax.set_xlabel('Satellite amount')
        ax.set_ylabel('Error [km]')
        ax.set_title(f'Error w.r.t {sat_amount} satellites for Newton-Gauss')
        bp = ax.boxplot(all_errors)  
        plt.setp(bp['whiskers'], color='k', linestyle='-')
        plt.setp(bp['fliers'], markersize=3.0)
        fig.savefig(f'figures/sat_num_{sat_amount}.png')

        if sat_amount == 5:
            continue

        max_errors.append(max_error)
        min_errors.append(min_error)
        mean_errors.append(mean_error)
        all_all_errors.append(all_errors)

    fig, ax = plt.subplots()
    ax.set_xlabel('Satellite amount')
    ax.set_ylabel('Error [km]')
    ax.set_title('Error w.r.t the amount of satellites for Newton-Gauss')
    bp = ax.boxplot(all_all_errors, positions=[i for i in range(6, sat_amount+1)])  
    plt.setp(bp['whiskers'], color='k', linestyle='-')
    plt.setp(bp['fliers'], markersize=3.0)
    fig.savefig(f'figures/all_sats.png')

    print([i for i in range(6, sat_amount+1)])
    print(max_errors)

    plt.clf()
    plt.plot([i for i in range(6, sat_amount+1)], max_errors, label="max")
    plt.plot([i for i in range(6, sat_amount+1)], min_errors, label="min")
    plt.plot([i for i in range(6, sat_amount+1)], mean_errors, label="mean")
    plt.xlabel("Satellite amount")
    plt.ylabel("Error [km]")
    plt.title("Different errors w.r.t. satellite amount with Newton-Gauss")
    plt.legend(loc='best')
    plt.xticks(np.arange(6, sat_amount+1, 1))
    plt.show()

main()