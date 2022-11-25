import numpy as np 
from Positions import Positions
import matplotlib.pyplot as plt

if __name__ == "__main__":

    x0 = np.array([0,0,6730,0])
    pos = Positions(x0, toll=1e-8)

    all_errors = []
    rand_thetas, rand_phis = pos.random_sat_angles(100,4)
    for i in range(0, len(rand_phis)):
        max_error = pos.distance_w_error(rand_thetas[i], rand_phis[i], err= 1e-8, tol_err=1e-8, sat_amount=4, method_type="Newton")
        all_errors.append(max_error)
        
    print('Max error in distance:', max(all_errors))
    print('Min error in distance:', min(all_errors))
    print('Average:', np.mean(all_errors))
    print('Standard deviation:', np.std(all_errors))
    fig = plt.figure(figsize =(10, 7))    
    plt.boxplot(all_errors)  
    plt.show()



