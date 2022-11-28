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
corr_phi = [pi/8, pi/6, 3*pi/8, pi/4]
corr_theta = [-pi/4, pi/2, 2*pi/3, pi/6]
incorr_phi = [pi/8 + 1e-8, pi/6 + 1e-8, 3*pi/8 - 1e-8, pi/4 - 1e-8]
A, B, C, t = [], [], [], []
for i in range(len(corr_phi)):
    t.append(location(corr_phi[i], corr_theta[i])['t']) #Vector of time for each sat t[s] derived from correct values
    values = location(incorr_phi[i], corr_theta[i]) #Derived from prerceived values
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
    eq_list = [0,0,0,0]
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

def main():
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    x,y,z,d = newtonmult(x0, 0.00001)
    print(f"Values: x = {x} ,  y = {y} , z = {z} , d =  {d}")
    print('The error is: ', end="")
    print(sqrt(pow(x - X, 2) + pow(y - Y, 2) + pow(z - Z, 2)), 'km')
    return 0
    

if __name__ == "__main__":
    main()