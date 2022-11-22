import matplotlib.pyplot as pst
import numpy as np 
from numpy import linalg as LA


# Constants
c = 299792.458 #Speed of light[km/s]
A = [15600, 18760, 17610, 19170]         #Vector of distances in plane A[km]
B = [7540, 2750, 14630, 610]             #Vector of distances in plane B[km]
C = [20140, 18610, 13480, 18390]         #Vector of distances in plane C[km]
t = [0.07074, 0.07220, 0.07690, 0.07242] #Vector of time for each sat t[s]


def F(x):
    '''Takes in a 4x1 vector as input for initial value. Creates a list and appends each equation 
    with relevant variable data and returns'''
    funcs = []
    for i in range(4):
        funcs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(c,2) * pow((t[i]-x[3]), 2))
    return funcs

def DF(x):
    """Creates each row of jacobi matrix independently(Hard coded) takes in a 4x1 vector
    of initial conditions and creates a matrix which is returned"""
    l1 = [(2*(x[0]-A[0])), (2*(x[1]-B[0])), (2*(x[2]-C[0])), (2*pow(c,2)* (t[0]-x[3]))]
    l2 = [(2*(x[0]-A[1])), (2*(x[1]-B[1])), (2*(x[2]-C[1])), (2*pow(c,2)* (t[1]-x[3]))]
    l3 = [(2*(x[0]-A[2])), (2*(x[1]-B[2])), (2*(x[2]-C[2])), (2*pow(c,2)* (t[2]-x[3]))]
    l4 = [(2*(x[0]-A[3])), (2*(x[1]-B[3])), (2*(x[2]-C[3])), (2*pow(c,2)* (t[3]-x[3]))]
    return np.array([l1,l2,l3,l4])

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

def main():
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    x,y,z,d = newtonmult(x0, 0.00000001)
    print("Final answers: \n")
    print("x is {:.6f} \t".format(x))
    print("y is {:.6f} \t".format(y))
    print("z is {:.6f} \t".format(z))
    print("d is {:.6f} \n".format(d))
    return 0 
main()