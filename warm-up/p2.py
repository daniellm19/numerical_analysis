import numpy as np

# constants
T = 1
n = 100
h = T/n
y_a = 1
y_b = -1

def nlbvpfd():
    w = np.zeros(n)
    for i in range(20):
        w = w - np.linalg.solve(jac(w), f(w))
    
    return w

def f(w):
    y = np.zeros(n)
    y[0] = y_a - w[0]*pow(h,2) + w[1]
    y[n-1] = w[n-2] - w[n-1]*pow(h, 2) + y_b
   
    
def jac(w):
    a = np.zeros((n,n))
    for i in range(1, n):
        a[i-1][i-1] = w[i-2] - w[i-1]*(2+pow(h, 2)) + w[i]
        
    for i in range(1, n):
        a[i-1][i] = 1
        a[i][i-1] = 1

w = nlbvpfd()
print(w)
    