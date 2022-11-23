import random
from numpy import pi


def random_angles():
    rand_phis = []
    rand_thetas = []
    for i in range(0,200):
        rand_phi = []
        rand_theta = []
        for j in range(0,4):
            rand_phi.append(random.uniform(0.0, pi/2))
            rand_theta.append(random.uniform(0.0, 2*pi))
        rand_thetas.append(rand_theta)
        rand_phis.append(rand_phi)
    return rand_thetas, rand_phis
