import numpy as np
from numpy import log, linalg as LA
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt


#constants
H = 0.005
K = 1.68
P = 5
L = 2
delta = 0.1

def mesh(xvals,yvals,w,xlabel='',ylabel='',zlabel=''):
    x,y = np.meshgrid(xvals,yvals)			# 3-D plot of 2D array w
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x, y, w, rstride=1)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_zlabel(zlabel)
    surf = ax.plot_surface(xvals,yvals,w, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def power_constant():
    return - P/(delta * K *L)

def const_1(k_or_h):
    return (-3/(2*k_or_h) + H/K)

def const_2(k_or_h):
    return (2/k_or_h)

def const_3(k_or_h):
    return -(1/(2*k_or_h))

def const_4(k,h):
    return -(2*((1/k)+(1/h)))

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
            print(y_loc)
            A[y_loc, (i-1)+j*m,] = const_5(k)      # The left node
            A[y_loc, (i+1)+j*m] = const_5(k)      # The right node 
            A[y_loc, i+j*m] = const_4(k,h)  # The node itself
            A[y_loc, i+(j-1)*m] = const_5(k)  # The bottom node
            A[y_loc, i+(j+1)*m] = const_5(k)  # The top node
            b[y_loc] = 0             # What the hell is this?
    plt.imshow(A)
    plt.colorbar()
    plt.show()
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
    plt.imshow(A)
    plt.colorbar()
    plt.show()
    # Þessi hluti er frekar sus, mundi definitely check'a á honum
    for j in range(n):	# left and right boundary points 
        i = 0
        A[i+j*m,i+j*m] = const_6(h)
        A[i+j*m, i+1+j*m] = -const_3(h)
        A[i+j*m, i+2+j*m] = const_3(h)
        b[i+j*m] = power_constant() # This assumes that all of the left side is power
        i = m-1
        A[i+j*m,i+j*m] = const_1(h)
        A[i+j*m, i-1+j*m] = const_2(h)
        A[i+j*m, i-2+j*m] = const_3(h)
        b[i+j*m] = 0 
    v = LA.solve(A,b)	# solve for solution in v labeling 
    w = np.reshape(v,(m,n),order='F') #translate from v to w

    print(f"A: {A}")
    print(f"b: {b}")

    plt.imshow(A)
    plt.colorbar()
    plt.show()
    mesh(x,y,w.T,'x','y','w') 
    return "ok"

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

w = poisson(0.2,0.1,0.2,0.1,4,4)
#A = fill_equation(5,5,1,2)
#print(A)