import random
from numpy import pi
import numpy as np
import json


def random_angles(pend_pairs_amount: int):
    '''Generates {pos} random positions for {sat_amount} satellites each'''
    rand_ini_vals = []
    for _ in range(pend_pairs_amount):
        y_0 = [random.uniform(pi/12, 2*pi), 0, random.uniform(pi/12, 2*pi), 0]
        rand_ini_vals.append(y_0)

    return rand_ini_vals

y_0 = random_angles(100)

with open("random_ini_values.json", "w") as f:
    f.write(json.dumps(y_0))


