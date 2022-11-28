import numpy as np 
from Positions import Positions

if __name__ == "__main__":
    '''Runs the program and gives stores the intitial guess
    And prints the solution in an acceptable way'''
    toll_error = 1e-6
    tiny_offset = 1e-5
    small_offset = 0.01

    tiny_phi = [np.pi/8-tiny_offset-tiny_offset, np.pi/8-tiny_offset, np.pi/8+tiny_offset, np.pi/8+tiny_offset+tiny_offset]
    tiny_theta = [np.pi/2-tiny_offset-tiny_offset, np.pi/2-tiny_offset, np.pi/2+tiny_offset, np.pi/2+tiny_offset+tiny_offset]

    small_phi = [np.pi/8-small_offset-small_offset, np.pi/8-small_offset, np.pi/8+small_offset, np.pi/8+small_offset+small_offset]
    small_theta = [np.pi/2-small_offset-small_offset, np.pi/2-small_offset, np.pi/2+small_offset, np.pi/2+small_offset+small_offset]
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method

    pos = Positions(x0, toll=toll_error)

    for i in range(0, len(tiny_phi)):
        pos.add_satellite_polar(tiny_phi[i], tiny_theta[i])

    x,y,z,d = pos.find_intersection()
    print(f"\nposition with e = {tiny_offset}: \n")
    print("x is {:.16f} \t".format(x))
    print("y is {:.16f} \t".format(y))
    print("z is {:.16f} \t".format(z))
    print("d is {:.16f} \n".format(d))

    pos.remove_all_sat()

    for i in range(0, len(small_phi)):
        pos.add_satellite_polar(small_phi[i], small_theta[i])

    x,y,z,d = pos.find_intersection()
    print(f"\nposition with e = {small_offset}: \n")
    print("x is {:.16f} \t".format(x))
    print("y is {:.16f} \t".format(y))
    print("z is {:.16f} \t".format(z))
    print("d is {:.16f} \n".format(d))