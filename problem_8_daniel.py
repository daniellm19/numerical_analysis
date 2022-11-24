import numpy as np
from numpy import sin, cos, sqrt, pi, mean, std, linalg as LA
import matplotlib.pyplot as plt
import random
from sympy.utilities.iterables import multiset_permutations
from operator import xor

# Constants
RHO = 26570                 # distance of satellites from earth[km]
LIGHT_SPEED = 299792.458    # speed of light [km/s]
X_0 = np.array([0,0,6370,0])          # Initial x, y, z, and d (time from sat to recv) values

def random_sat_angles(pos: int, rows: int):
    '''
    Generates {pos} random positions for {rows} amount of satellites each
    
    Arguments
    ----------
    pos : The number of random positions to calculate for the satellites
    rows : The amount of satellites that become rows in f(x) and df(x) (Jacobi)
    '''
    rand_phis = []
    rand_thetas = []
    for _ in range(0, pos):
        rand_phi = []
        rand_theta = []
        for _ in range(0, rows):
            rand_phi.append(random.uniform(0.0, pi/2))
            rand_theta.append(random.uniform(0.0, 2*pi))
        rand_thetas.append(rand_theta)
        rand_phis.append(rand_phi)
    return rand_thetas, rand_phis

def cartesian_calc(phi: float, theta: float) -> dict:
    """
    Calculates the cartesian coordinates of the satellites along with the time it takes
    to reach the satellite from the inital position of the receiver

    Arguments
    ----------
    phi : The altitude of the satellite (between 0 and pi/2 for simplification)
    theta : The polar angle of the satellite (between 0 and 2*pi)
    """
    if (phi>(np.pi/2)) or (0>phi):
        raise ValueError("Phi (the altitude) should be between 0 and pi/2")
    if (theta>(2*np.pi)) or (0>theta):
        raise ValueError("Theta (the polar angle) should be between 0 and 2*pi")

    A = RHO * sin(phi) * cos(theta) if abs(RHO * sin(phi) * cos(theta)) > 1e-10 else 0
    B = RHO * sin(phi) * sin(theta) if abs(RHO * sin(phi) * sin(theta)) > 1e-10 else 0
    C = RHO * cos(phi) if abs(RHO * cos(phi)) > 1e-10 else 0
    distance = sqrt(pow((X_0[0] - A), 2) + pow((X_0[1] - B), 2) + pow((X_0[2] - C), 2))
    t = distance / LIGHT_SPEED

    ret_dict = {'A': A, 'B': B, 'C': C, 'distance': distance, 't': t}
    return ret_dict

def get_incorr_angle(err: int, angles: list) -> list:
    """
    Creates all of the combinations of a error value being either added or subtracted from
    a list of angle coordinates
    
    Arguments
    ----------
    err : The error that should be used for the permutation
    angles : The list of angles that should be error-riddled
    """
    all_err_perm = []
    err_list = [err] * len(angles)
    for i in range(0,len(err_list)+1):
        err_copy = err_list.copy()
        for j in range(0,i):
            err_copy[j] = -err_list[j]
        multiset_perm = multiset_permutations(err_copy)
        for i in multiset_perm:
            all_err_perm.append(i)

    # Adds all of the error permutations to the a_list
    all_perms = [[a + b for a, b in zip(all_err_perm[i], angles)] for i in range(len(all_err_perm))]

    return all_perms

def calculate_cartesian(theta: float, phi: float) -> float: #TODO indicate four float value return types
    """
    Gets the cartesian values and the distance from the angles (altitude and azimuth/polar) from
    the satellites

    Arguments
    ----------
    theta : The altitude angles of the satellites
    phi : The azimuth/polar angles of the satellites
    """
    if (len(theta) != len(phi)) or (len(theta)==0) or (len(phi)==0):
        raise ValueError("The thetas and phis should not be zero and should have the same length")

    A, B, C, t = [], [], [], []
    for i in range(len(phi)):
        values = cartesian_calc(phi[i], theta[i])     #Derived from prerceived values
        A.append(values['A'])                        #Vector of distances in plane A[km]
        B.append(values['B'])                        #Vector of distances in plane B[km]
        C.append(values['C'])                        #Vector of distances in plane C[km]
        t.append(values['t'])                        #Vector of time for each sat t[s] derived from correct values
    return A, B, C, t

def f(x: list, A: list, B: list, C: list, t: list, rows: int) -> list: #TODO have return type be np.array
    '''
    Takes in a 1x4 vector of x, y, z, and d (t_{recv}-t_{sat}) for the initial value along with the
    cartesian positions of the satellites. 
    Creates a matrix and appends each equation with relevant variable data into the list 
    and returns a {rows}x1 matrix of the equations.

    Arguments
    ----------
    x : A 1x4 vector of x, y, z, and d (t_{recv}-t_{sat})
    A : A list of the x-axis coordinates of the satellites
    B : B list of the x-axis coordinates of the satellites
    C : C list of the x-axis coordinates of the satellites
    t : The time it takes to travel from the satellite to the receiver (at light speed)
    rows : The amount of satellites (which translates to equations)
    '''
    funcs = []
    for i in range(rows):
        funcs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(LIGHT_SPEED,2) * pow((t[i]-x[3]), 2))
    return funcs

def df(x: list, A: list, B: list, C: list, t: list, rows: int) -> list: #TODO have return type be np.array
    '''
    Calculates each row of jacobi matrix independently (Hard coded). Takes in a 5x1 vector
    of initial conditions and creates a {rows}x{rows} matrix
    
    Arguments
    ----------
    x : A 1x4 vector of x, y, z, and d (t_{recv}-t_{sat})
    A : A list of the x-axis coordinates of the satellites
    B : B list of the x-axis coordinates of the satellites
    C : C list of the x-axis coordinates of the satellites
    t : The time it takes to travel from the satellite to the receiver (at light speed)
    rows : The amount of satellites (which translates to equations)
    '''
    eq_list = [0]*rows
    for i in range(0, len(eq_list)):
        eq_list[i] = [(2*(x[0]-A[i])), (2*(x[1]-B[i])), (2*(x[2]-C[i])), (2*pow(LIGHT_SPEED,2)* (t[i]-x[3]))]
    return np.array(eq_list)

def newton_mult(x0: list, tol: int, A: list, B: list, C: list, t: list, rows: int) -> list:
    '''
    Calculates the intercept that the satellites estimate with the Newton method.
    Takes in inital values along with a tolerance to determine when the root has
    been reached.

    Arguments
    ----------
    x0 : The inital estimate of where the intercept root is located
    tol : The tolerance of when the root has been reached
    A : A list of the x-axis coordinates of the satellites
    B : B list of the x-axis coordinates of the satellites
    C : C list of the x-axis coordinates of the satellites
    t : The time it takes to travel from the satellite to the receiver (at light speed)
    rows : The amount of rows (satellites)
    '''
    x=x0
    oldx=x+2*tol
    while LA.norm(x-oldx, np.inf)>tol:
        oldx=x
        s=-LA.solve(df(x, A, B, C, t, rows), f(x, A, B, C, t, rows))
        x=x+s
    return(x)

def newton_gauss_mult(x0: list, tol: int, A: list, B: list, C: list, t: list, rows: int) -> list:
    """
    Calculates and returns the estimated intercept of the receivers position using the Newton-Gauss method.

    Arguments
    ----------
    x0 : The inital estimate of where the intercept root is located
    tol : The tolerance of when the root has been reached
    A : A list of the x-axis coordinates of the satellites
    B : B list of the x-axis coordinates of the satellites
    C : C list of the x-axis coordinates of the satellites
    t : The time it takes to travel from the satellite to the receiver (at light speed)
    rows : The amount of rows (satellites)
    """
    x = x0
    oldx=x+2*tol
    while LA.norm(x-oldx, np.inf)>tol:
        oldx = x
        left_side = np.dot(np.transpose(df(x, A, B, C, t, rows)), df(x, A, B, C, t, rows))
        right_side = np.dot((np.transpose(df(x, A, B, C, t, rows))), f(x, A, B, C, t, rows))
        s = -LA.solve(left_side, right_side)
        x = x+s
    return(x)

def distance_w_error(thetas: list, phis: list, err: int, tol_err: int, sat_amount: int, method_type: str = "NG") -> float:
    '''
    Takes a list of angle coordinates of satellites and adds errors to those coordinates,
    then it calculates the error (using either Newton or Newton-Gauss) and returns 
    the maximum error found

    Arguments
    ----------
    thetas : The altitude angles
    phis: The azimuth/polar angles
    err : The error you want to add to the satellites
    tol_err : The tolerance error for either the Newton or Newton-Gauss method
    sat_amount : The amount of satellites
    method_type: Specifies either the Newton method ("G") or the Newton-Gauss method ("NG")
    t : The time it takes to travel from the satellite to the receiver (at light speed)
    rows : The amount of rows (satellites)
    '''

    if (method_type!="NG") and (method_type!="N"):
        raise ValueError("method_type should be either 'NG' (for Newton-Gauss) or 'N' (for Newton)")

    incorr_phis = get_incorr_angle(err, phis)
    all_errors = []
    for incorr_phi in incorr_phis:
        A, B, C, _ = calculate_cartesian(thetas, incorr_phi)
        _, _, _, t = calculate_cartesian(thetas, phis)
        if method_type=="NG":
            x,y,z,_ = newton_gauss_mult(X_0, tol_err, A, B, C, t, sat_amount)
        elif method_type=="N":
            x,y,z,_ = newton_mult(X_0, tol_err, A, B, C, t, sat_amount)
        all_errors.append(sqrt(pow(x - X_0[0], 2) + pow(y - X_0[1], 2) + pow(z - X_0[2], 2)))
    return max(all_errors)

def main():
    sat_pos_amount = 100 
    min_sat_amount = 5
    sat_amount = 9
    all_all_errors = []

    for sat_amount in range(min_sat_amount, sat_amount+1):
        all_errors = []

        rand_thetas, rand_phis = random_sat_angles(sat_pos_amount, sat_amount)
        for i in range(len(rand_phis)):
            max_error = distance_w_error(rand_thetas[i], rand_phis[i], 1e-8, 1e-8, sat_amount, method_type="NG")
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

        if sat_amount <= min_sat_amount:
            continue
        all_all_errors.append(all_errors)

    fig, ax = plt.subplots()
    ax.set_xlabel('Satellite amount')
    ax.set_ylabel('Error [km]')
    ax.set_title('Error w.r.t the amount of satellites for Newton-Gauss')
    bp = ax.boxplot(all_all_errors, positions=[i for i in range(6, sat_amount+1)])  
    plt.setp(bp['whiskers'], color='k', linestyle='-')
    plt.setp(bp['fliers'], markersize=3.0)
    fig.savefig(f'figures/all_sats.png')

    plt.clf()
    plt.plot([i for i in range(6, sat_amount+1)], [np.max(i) for i in all_all_errors], label="max")
    plt.plot([i for i in range(6, sat_amount+1)], [np.min(i) for i in all_all_errors], label="min")
    plt.plot([i for i in range(6, sat_amount+1)], [np.mean(i) for i in all_all_errors], label="mean")
    plt.xlabel("Satellite amount")
    plt.ylabel("Error [km]")
    plt.title("Different errors w.r.t. satellite amount with Newton-Gauss")
    plt.legend(loc='best')
    plt.xticks(np.arange(6, sat_amount+1, 1))
    plt.savefig('figures/diff_errors')
    plt.show()

main()