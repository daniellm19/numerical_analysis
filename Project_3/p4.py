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

def power_constant():
    return  -P/(delta * K *L)

def const_1(k_or_h):
    return (-3/(2*k_or_h) + H/K)

def const_2(k_or_h):
    return (2/k_or_h)

def const_3(k_or_h):
    return -(1/(2*k_or_h))

def const_4(k,h):
    return (-2*((1/pow(k,2))+(1/pow(h,2))+(H/(K*delta))))

def const_5(k_or_h):
    return (1/pow(k_or_h,2))

def const_6(h):
    return (-3/(2*h))

def poisson(xl,xr,yb,yt,M,N):
    m = M
    n = N
    h = (xr-xl)/M
    k = (yt-yb)/N
    x = np.linspace(xl,xr,m)    # set mesh values 
    y = np.linspace(yb,yt,n)
    A = np.zeros((m*n,m*n))     # The A matrix is the reverse of what is seen in the
                                # lectures because numpy's indexing is insane and
                                # I won't bother to change it
    b = np.zeros(m*n) 
    for i in range(1,m-1):  # interior points 
        for j in range(1,n-1):
            y_loc = i+j*m
            A[y_loc, (i-1)+j*m,] = const_5(h)      # The left node
            A[y_loc, (i+1)+j*m] = const_5(h)      # The right node 
            A[y_loc, i+j*m] = const_4(k,h)  # The node itself
            A[y_loc, i+(j-1)*m] = const_5(k)  # The bottom node
            A[y_loc, i+(j+1)*m] = const_5(k)  # The top node
            b[y_loc] = 0             # What the hell is this?
    #plt.imshow(A)
    #plt.colorbar()
    #plt.show()
    for i in range(1,m-1): 		# bottom and top boundary points 
        # You should look at i like 1 here and j like 1 as well
        j = n-1               # Top
        A[i+j*m, i+j*m] = const_1(k)  # Why is it set to one?
        A[i+j*m, i+(j-1)*m] = const_2(k)  # This needs to be error checked
        A[i+j*m, i+(j-2)*m] = const_3(k)  # This needs to be error checked
        b[i+j*m] = 0 # I think this doesn't make sense, we already have all boundary conditions
        j = 0             # Bottom
        A[i+j*m,i+j*m] = const_1(k)  # Why?
        A[i+j*m, i+(j+1)*m] = const_2(k)  # This needs to be error checked
        A[i+j*m, i+(j+2)*m] = const_3(k)  # This needs to be error checked
        b[i+j*m] = 0 # Why?
    #plt.imshow(A)
    #plt.colorbar()
    #plt.show()
    # Þessi hluti er frekar sus, mundi definitely check'a á honum
    for j in range(n):	# left and right boundary points 
        i = 0
        A[i+j*m,i+j*m] = const_6(h)
        A[i+j*m, i+1+j*m] = const_2(h)
        A[i+j*m, i+2+j*m] = const_3(h)
        b[i+j*m] = power_constant() # This assumes that all of the left side is power
        i = m-1
        A[i+j*m,i+j*m] = const_1(h)
        A[i+j*m, i-1+j*m] = const_2(h)
        A[i+j*m, i-2+j*m] = const_3(h)
        b[i+j*m] = 0 
    #plt.imshow(A)
    #plt.colorbar()
    #plt.show()
    v = LA.solve(A,b)	# solve for solution in v labeling 
    v = [i-20 for i in v]
    w = np.reshape(v,(m,n),order='F') #translate from v to w
    return w

def fill_equation(n,m,low,high):
    A = np.zeros((n,m))
    for i in range(0,n):    # The y-axis, 0 is the top, n-1 is the bottom
        for j in range(0,m):# The x-axis, 0 is the left, n-1 is the right
            if i==0:
                A[i,j] = 3  # Top side
                continue
            elif i==(m-1):
                A[i,j] = 4  # Bottom side
                continue
            elif j==0:
                if ((low<=i) and (i<=high)):
                    A[i,j] = 6  # Power side
                    continue
                else:
                    A[i,j] = 1  # Left side
                    continue
            elif j == (n-1):
                A[i,j] = 2  # Right side
                continue
            else:
                A[i,j] = 5  # Inside
    
    return A

def deviation(reference, w):
    print('hæ')
    return 1

def main():
    ns = ms = list(range(20, 100, 10))
    reference_w = poisson(0.2,0.1,0.2,0.1,100,100)
    too_long, too_much_error, good_values = [], [], []
    for n in ns:
        print
        for m in ms:
            start_time = time()
            w = poisson(0.2,0.1,0.2,0.1,m,n)
            duration = time() - start_time
            print(duration)
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