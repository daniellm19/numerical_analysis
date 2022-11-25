import numpy as np 
from numpy import sin, cos, sqrt, pi, linalg as LA
import matplotlib.pyplot as pst

rho = 26570 # kilometers
c = 299792.458 # speed of light [km/s]
X = 0
Y = 0
Z = 6370

def location(phi, theta):
    A = rho * sin(phi) * cos(theta) if abs(rho * sin(phi) * cos(theta)) > 1e-10 else 0
    B = rho * sin(phi) * sin(theta) if abs(rho * sin(phi) * sin(theta)) > 1e-10 else 0
    C = rho * cos(phi) if abs(rho * cos(phi)) > 1e-10 else 0
    distance = sqrt(pow((X - A), 2) + pow((Y - B), 2) + pow((Z - C), 2))
    t = distance / c

    ret_dict = {'A':A, 'B':B, 'C':C, 'distance':distance, 't':t}
    #print('A =', A, 'km\nB =', B, 'km\nC =', C, 'km\nDistance from north pole =', distance, 'km\nt =', t, 'sec\n\n')
    return ret_dict

#Problem 3 starts here:

#Constants
tol_error = 1e-6
tiny_offset = 1e-5
small_offset = 0.01


tiny_phi = [pi/8-tiny_offset-tiny_offset, pi/8-tiny_offset, pi/8+tiny_offset, pi/8+tiny_offset+tiny_offset]
tiny_theta = [pi/2-tiny_offset-tiny_offset, pi/2-tiny_offset, pi/2+tiny_offset, pi/2+tiny_offset+tiny_offset]

small_phi = [pi/8-small_offset-small_offset, pi/8-small_offset, pi/8+small_offset, pi/8+small_offset+small_offset]
small_theta = [pi/2-small_offset-small_offset, pi/2-small_offset, pi/2+small_offset, pi/2+small_offset+small_offset]

A, B, C, t = [], [], [], []
def populate_matrix(phi, theta):
    A.clear(), B.clear(), C.clear(), t.clear()
    for i in range(len(phi)):
        t.append(location(phi[i], theta[i])['t']) #Vector of time for each sat t[s] derived from correct values
        values = location(phi[i], theta[i]) #Derived from prerceived values
        A.append(values['A']) #Vector of distances in plane A[km]
        B.append(values['B']) #Vector of distances in plane B[km]
        C.append(values['C']) #Vector of distances in plane C[km]

def F(x):
    funcs = []
    for i in range(4):
        funcs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(c,2) * pow((t[i]-x[3]), 2))
    return funcs

def DF(x):
    '''Creates each row of jacobi matrix independently(Hard coded) takes in a 4x1 vector
    of initial conditions and creates a matrix which is returned'''
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
    populate_matrix(tiny_phi, tiny_theta)
    x,y,z,d = newtonmult(x0, tol_error)
    print(f"\nposition with e = {tiny_offset}: \n")
    print("x is {:.16f} \t".format(x))
    print("y is {:.16f} \t".format(y))
    print("z is {:.16f} \t".format(z))
    print("d is {:.16f} \n".format(d))
    populate_matrix(small_phi, small_theta)
    x,y,z,d = newtonmult(x0, tol_error)
    print(f"\nposition with e = {small_offset}: \n")
    print("x is {:.16f} \t".format(x))
    print("y is {:.16f} \t".format(y))
    print("z is {:.16f} \t".format(z))
    print("d is {:.16f} \n".format(d))
    return 0 
main()

if __name__ == "__main__":
    print(main)