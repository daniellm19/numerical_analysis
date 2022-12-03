import matplotlib.pyplot as plt
import numpy as np

def x(t, x_0, w):
    return x_0*np.cos(w*t)
def x_prime(t, x_0, w):
    return -x_0*w*np.sin(w*t)

t_range = np.arange(0, 2*np.pi, np.pi/4)
for t in t_range:
    plt.plot(x(t, 1, 1), x_prime(t, 1, 1), markersize=3, marker='o')

a = x_prime(t_range, 1, 1)
b = x(t_range, 1, 1)
a = np.append(a, [a[0]])
b = np.append(b, [b[0]])

plt.plot(a, b)
plt.show()