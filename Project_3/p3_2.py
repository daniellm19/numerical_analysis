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

def poisson(xl,xr,yb,yt,M,N):
    mn=(M+1)*(N+1)
    m = M+1
    n = N+1
    h = L/(M)
    k = L/(N)
    print(f"h: {h}")
    print(f"k: {k}")
    x = np.linspace(xl,xr,m)    # set mesh values 
    y = np.linspace(yb,yt,n)
    A = np.zeros((mn,mn))
    h2 = pow(h,2)
    k2 = pow(k,2)
    b = np.zeros((mn,1)) 

    print(f"1/h2: {1/h2}")
    print(f"1/k2: {1/k2}")
    print(f"-((2/k2)+(2/h2)+((2*H)/(K*delta))): {-((2/k2)+(2/h2)+((2*H)/(K*delta)))}")
    print(f"((-3/(2*k)) + (H/K)): {((-3/(2*k)) + (H/K))}")
    print(f"(2/k): {(2/k)}")
    print(f"(-P)/(delta*K*L): {(-P)/(delta*K*L)}")

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
        A[y_loc, y_loc] =       ((-3/(2*k)) + (H/K))
        A[y_loc, y_loc+m] =     (2/k)  
        A[y_loc, y_loc+m*2] =   (-1/(2*k))  

        j = n               # Top
        y_loc = i+(j-1)*m-1
        A[y_loc, y_loc] =       ((-3/(2*k)) + (H/K))
        A[y_loc, y_loc-m] =     (2/k)
        A[y_loc, y_loc-m*2] =   (-1/(2*k))  

    for j in range(1,n+1):	# left and Right boundary points 
        i = 1               # Power (left)
        y_loc = i+(j-1)*m-1
        A[y_loc, y_loc] =   -3
        A[y_loc, y_loc+1] = 4
        A[y_loc, y_loc+2] = -1
        b[y_loc] = (-2*h*P)/(delta*K*L)

        i = m               # Right
        y_loc = i+(j-1)*m-1
        A[y_loc, y_loc] =   (-3 + ((2*h*H)/K))
        A[y_loc, y_loc-1] = 4
        A[y_loc, y_loc-2] = -1

    v = LA.solve(A,b)	    # solve for solution in v labeling 
    v = [i+20 for i in v]
    w = np.reshape(v,(m,n),order='F') #translate from v to w

    print(f"A: {A}")
    print(f"b: {b}")

    print(f"v: {v}")
    print(f"max v: {max(v)}")
    print(f"w: {w}")

    fig, ax = plt.subplots()

    c = ax.pcolormesh(x, y, w, cmap='RdBu', vmin=w.min(), vmax=w.max())
    ax.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)

    plt.show()

    plt.imshow(A)
    plt.colorbar()
    plt.show()
    mesh(x,y,w.T,'x','y','w') 
    return "ok"

w = poisson(0,0.2,0,0.2,9,9)
#w = poisson(0,1,1,2,4,4)