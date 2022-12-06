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

def f(x,y):
    return 0

def g1(x):
    return -(H/K) * x

def g2(x): 
    return -(H/K) * x

def g3(y): # Left
    return - P/(delta * K *L)

def g4(y): 
    return - (H/K) * y


def poisson(xl,xr,yb,yt,M,N):
    m = M+1
    n = N+1 
    h = (xr-xl)/M
    k = (yt-yb)/N
    x = np.linspace(xl,xr,m) # set mesh values 
    y = np.linspace(yb,yt,n)
    A = np.zeros((m*n,m*n))
    b = np.zeros(m*n) 
    for i in range(1,m-1):  # interior points 
        for j in range(1,n-1):
            A[i+j*m,i-1+j*m] = 1/pow(h, 2)
            A[i+j*m,i+1+j*m] = 1/pow(h, 2)
            A[i+j*m,i  +j*m] = (2*H)/(K*delta)
            A[i+j*m,i  +(j-1)*m] = 1/pow(k, 2)
            A[i+j*m,i  +(j+1)*m] = 1/pow(k, 2)
            b[i+j*m] = f(x[i],y[j]) 
    for i in range(m): 		# bottom and top boundary points 
        j = 0  
        A[i+j*m,i+j*m] = 1
        b[i+j*m] = g1(x[i]) 
        j = n-1
        A[i+j*m,i+j*m] = 1
        b[i+j*m] = g2(x[i]) 
    for j in range(1,n-1):	# left and right boundary points 
        i = 0
        A[i+j*m,i+j*m] = 1
        b[i+j*m] = g3(y[j]) 
        i = m-1
        A[i+j*m,i+j*m] = 1
        b[i+j*m] = g4(y[j]) 
    v = LA.solve(A,b)	# solve for solution in v labeling 
    w = np.reshape(v,(m,n),order='F') #translate from v to w
    print(w)
    mesh(x,y,w.T,'x','y','w') 
    return w

poisson(0,2,0,2,10,10) 