import math
import matplotlib.pyplot as pst
import numpy as np 
from numpy import sin, cos, sqrt, pi, linalg as LA

rho = 26570
error = 10e-8

phis = [pi/8, pi/6, 3*pi/8, pi/4]
thetas = [-pi/4, pi/2, 2*pi/3, pi/6]

error_phis = [pi/8+error, pi/6+error, 3*pi/8-error, pi/4-error]
error_thetas = thetas

c = 299792.458 #Speed of light[km/s]
A = []
B = []
C = []
t = []

def calculate_abc(phis: list, thetas: list):
    for i in range(0,len(phis)):
        A_val = rho*sin(phis[i])*cos(thetas[i])
        A.append(A_val)
        B_val = rho*sin(phis[i])*sin(thetas[i])
        B.append(B_val)
        C_val = rho*cos(phis[i])
        C.append(C_val)

def calculate_t(phis: list, thetas: list):
    for i in range(0,len(phis)):
        A_val = rho*sin(phis[i])*cos(thetas[i])
        B_val = rho*sin(phis[i])*sin(thetas[i])
        C_val = rho*cos(phis[i])
        t_val = sqrt(pow(A_val,2)+pow(B_val,2)+pow(C_val,2))/C_val
        t.append(t_val)


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
    calculate_abc(error_phis, error_thetas)
    calculate_t(phis, thetas)
    print(A)
    print(B)
    print(C)
    print(t)
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    x,y,z,d = newtonmult(x0, 0.1)
    print("Final answers: \n")
    print("x is {:.6f} \t".format(x))
    print("y is {:.6f} \t".format(y))
    print("z is {:.6f} \t".format(z))
    print("d is {:.6f} \n".format(d))
    return 0 
main()