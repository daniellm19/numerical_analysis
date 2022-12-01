import numpy as np
import matplotlib.pyplot as plt


ax = plt.figure().add_subplot()

# Prepare arrays x, y, z
theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
x = np.linspace(-2, 2, 100)
r = x**2 + 1
y = r * np.sin(theta)

ax.plot(x, y, label='parametric curve')
ax.legend()

plt.show()