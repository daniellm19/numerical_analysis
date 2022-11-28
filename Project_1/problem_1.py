import numpy as np 
from Positions import Positions

def main():
    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    pos = Positions(pos_guess=x0, toll=1e-8)
    A = [15600, 18760, 17610, 19170]         # Vector of distances in plane A[km]
    B = [7540, 2750, 14630, 610]             # Vector of distances in plane B[km]
    C = [20140, 18610, 13480, 18390]         # Vector of distances in plane C[km]
    t = [0.07074, 0.07220, 0.07690, 0.07242] # Vector of time for each sat t[s]
    for i in range(0,4):
        pos.add_satellite_cart(A[i], B[i], C[i], t[i])
    x,y,z,d = pos.find_intersection()
    print("Final answers: \n")
    print("x is {:.4f} \t".format(x))
    print("y is {:.4f} \t".format(y))
    print("z is {:.2f} \t".format(z))
    print("d is {:.8f} \n".format(d))
    return 0 
main()