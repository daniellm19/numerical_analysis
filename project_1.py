import numpy as np
from numpy import linalg as LA

C = 299792.458

const_val = np.array([[15600, 7540, 20140, 0.07074],[18760, 2750, 18610, 0.07220],[17610, 14630, 13480, 0.07690],[19170, 610, 18390, 0.07242]])

def F(x):
    f_array = np.array([[0], [0], [0], [0]])

    for i in range(0,4):
        f_array[i][0] = pow((x[0]-const_val[i][0]),2)+pow((x[1]-const_val[i][1]),2)+pow((x[2]-const_val[i][2]),2)-pow(C,2)*pow((const_val[i][3]-x[3]),2)

    return f_array

def DF(x):
    f_array = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]

    for i in range(0,4):
        f_array[i][0] = 2*(x[0]-const_val[i][0])
        f_array[i][1] = 2*(x[1]-const_val[i][1])
        f_array[i][2] = 2*(x[2]-const_val[i][2])
        f_array[i][3] = -2*pow(C,2)*(const_val[i][3]-x[3])

    return np.array(f_array)

def newton_mult(x0,tol):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert ráð fyrir að F(x) og Jacobi fylki DF(x) séu skilgreind annars staðar'''
    x=x0
    oldx=x+2*tol
    while LA.norm(x-oldx,np.inf)>tol:
        oldx=x
        print(x)
        s=-LA.solve(DF(x),F(x))
        print(s)
        x=x+s
        print(x)
    return(x)
        
if "__main__" == __name__:
    newton_mult(np.array([0,0,6370,0]),0.0000001)