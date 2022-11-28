import numpy as np 
from Positions import Positions

if __name__ == "__main__":
    pos = Positions(np.array([0,0,6730,0]), 1e-8)
    print(pos.cartesian_calc(np.pi/3, np.pi/2))