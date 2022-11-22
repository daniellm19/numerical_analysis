import numpy as np
from numpy import linalg as LA

def F(x):
    return pow(x,3)-4*pow(x,2)+3*x+2

def DF(x):
    return 3*pow(x,2)-8*x+3

def newton_mult(x0,tol):
    '''x0 er vigur i R^n skilgreindur t.d. sem
    x0=np.array([1,2,3])
    gert ráð fyrir að F(x) og Jacobi fylki DF(x) séu skilgreind annars staðar'''
    x=x0
    oldx=x+2*tol
    while LA.norm(x-oldx,np.inf)>tol:
        oldx=x
        s=-LA.solve(DF(x),F(x))
        x=x+s
    return(x)
        
if "__main__" == __name__:
    newton_mult(np.array([1,2,3]),0.0000001)