import random
from numpy import pi
import json


def random_angles(pos: int, sat_amount: int):
    '''Generates {pos} random positions for {sat_amount} satellites each'''
    rand_phis = []
    rand_thetas = []
    for i in range(0, pos):
        rand_phi = []
        rand_theta = []
        for j in range(0, sat_amount):
            rand_phi.append(random.uniform(0.0, pi/2))
            rand_theta.append(random.uniform(0.0, 2*pi))
        rand_thetas.append(rand_theta)
        rand_phis.append(rand_phi)
    return rand_thetas, rand_phis

thetas, phis = random_angles(100, 4)
all_angles = [thetas, phis]

with open("random_angles.json", "w") as f:
    f.write(json.dumps(all_angles))


