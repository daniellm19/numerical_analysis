import numpy as np
import matplotlib.pyplot as plt


def nlbvpfd(inter, bv, n):
    [a,b] = inter; [ya,yb] = bv
    h = (b-a)/(n+1)
    w = np.zeros(n)				# Initial guess for Newton
    old_w = w + 2 * 1e-5
    while np.linalg.norm(w-old_w, np.inf)>1e-5:	# Newton
        old_w = w
        print(w)
        s = - np.linalg.solve(jac(w, h), f(w, h))
        w += s


    plt.plot( np.linspace(a,b,n+2), np.hstack((ya,w,yb)), linewidth=2)
    return w

def f(w, h):
    n = len(w)
    y = np.zeros(n)
    y[0] = 1
    y[-1] = -1
    for i in range(1,n-1):
        y[i] =  w[i-1] - w[i] * (2 + pow(h, 2)) + w[i+1]
    return y

def jac(w, h):
	n = len(w)
	a = np.zeros((n,n))
	for i in range(n):
		a[i,i] = - (2 + pow(h, 2))
	for i in range(n-1):
		a[i,i+1] = 1
		a[i+1,i] = 1
	return a

n = 40
a = 0
b = 1
w = nlbvpfd([a,b],[1,-1],  n)
print(w)

plt.show()