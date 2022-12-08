import numpy as np
from numpy import log, linalg as LA
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
from time import time

#constants
H = 0.005
K = 1.68
P = 5
L = 2
delta = 0.1

def poisson(xl,xr,yb,yt,M,N):
    mn=(M+1)*(N+1)
    m = M+1
    n = N+1
    Lx = (xr-xl)
    Ly = (yt-yb)
    h = Lx/(M)
    k = Ly/(N)
    x = np.linspace(xl,xr,m)    # set mesh values 
    y = np.linspace(yb,yt,n)
    A = np.zeros((mn,mn))
    h2 = pow(h,2)
    k2 = pow(k,2)
    b = np.zeros((mn,1)) 

    for i in range(2,m):    # interior points 
        for j in range(2,n):
            y_loc = i+(j-1)*m-1
            A[y_loc, y_loc-1] = 1/h2        # The left node
            A[y_loc, y_loc+1] = 1/h2        # The right node 
            A[y_loc, y_loc] =   -((2/k2)+(2/h2)+((2*H)/(K*delta)))  # The node itself
            A[y_loc, y_loc-m] = 1/k2        # The bottom node
            A[y_loc, y_loc+m] = 1/k2        # The top node      

    for i in range(2,m):    # Bottom and Top boundary points 
        j = 1               # Bottom
        y_loc = i+(j-1)*m-1
        A[y_loc, y_loc] =       (-3 + (2*h*H/K))
        A[y_loc, y_loc+m] =     4  
        A[y_loc, y_loc+m*2] =   -1

        j = n               # Top
        y_loc = i+(j-1)*m-1
        A[y_loc, y_loc] =       (-3 + ((2*h*H)/K))
        A[y_loc, y_loc-m] =     4
        A[y_loc, y_loc-m*2] =   -1  

    for j in range(1,n+1):	# left and Right boundary points 
        i = 1               # Power (left)
        y_loc = i+(j-1)*m-1
        A[y_loc, y_loc] =       -3
        A[y_loc, y_loc+1] =     4
        A[y_loc, y_loc+2] =     -1
        b[y_loc] = (-2*h*P)/(delta*K*L)

        i = m               # Right
        y_loc = i+(j-1)*m-1
        A[y_loc, y_loc] =       (-3 + ((2*h*H)/K))
        A[y_loc, y_loc-1] =     4
        A[y_loc, y_loc-2] =     -1

    v = LA.solve(A,b)	    # solve for solution in v labeling 
    v = [i+20 for i in v]
    w = np.reshape(v,(m,n),order='F') #translate from v to w
    return w


def deviation(reference, w):
    difference = [abs(w[0][0] - reference[0][0]), abs(w[-1][0] - reference[-1][0]), abs(w[0][-1] - reference[0][-1]), abs(w[-1][-1] - reference[-1][-1])]
    return max(difference)

def main():
    ns = ms = list(range(20, 100, 10))
    reference_w = poisson(0.2,0.1,0.2,0.1,100,100)
    too_long, too_much_error, good_values = [], [], []
    for n in ns:
        print(n)
        for m in ms:
            start_time = time()
            w = poisson(0.2,0.1,0.2,0.1,m,n)
            duration = time() - start_time
            max_deviation = deviation(reference_w, w)
            if duration >= 0.5:
                too_long.append({'m': m, 'n': n, 'Time': duration})
            if max_deviation >= 0.01:
                too_much_error.append({'m': m, 'n': n, 'Max error': max_deviation})
            if duration < 0.5 and max_deviation < 0.01:
                good_values.append({'m': m, 'n': n, 'Time': duration, 'Max error': max_deviation})
    
    print('Good values:', good_values)
    print('Takes too much time:', too_long)
    print('Too much error:', too_much_error)


main()