import numpy as np 
from Positions import Positions

if __name__ == "__main__":
    corr_phi = [np.pi/8, np.pi/6, 3*np.pi/8, np.pi/4]
    corr_theta = [-np.pi/4, np.pi/2, 2*np.pi/3, np.pi/6]

    x0 = np.array([0,0,6370,0]) #Initial guess for newtons method
    pos = Positions(x0, toll=1e-8)
    incorr_phis = pos.get_incorr_angle(err= 1e-8, angles = corr_phi)

    all_errors = []
    for incorr_phi in incorr_phis:
        for i in range(0, len(incorr_phi)):
            t = pos.cartesian_calc(corr_phi[i], corr_theta[i])['t']         #Vector of time for each sat t[s] derived from correct values
            abc_values = pos.cartesian_calc(incorr_phi[i], corr_theta[i])   #Derived from perrceived values
            pos.add_satellite_cart(abc_values['A'], abc_values['B'], abc_values['C'], t)

        x, y, z, _ = pos.find_intersection()
        all_errors.append(np.sqrt(pow(x - pos.x0[0], 2) + pow(y - pos.x0[1], 2) + pow(z - pos.x0[2], 2)))
        pos.remove_all_sat()
    print('The maximum error is:', max(all_errors))