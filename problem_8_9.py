import numpy as np
from numpy import sin, cos, sqrt, pi, mean, std, linalg as LA
import matplotlib.pyplot as plt
import random
from sympy.utilities.iterables import multiset_permutations
from operator import xor
from Positions import Positions


if __name__ == "__main__":
    sat_pos_amount = 100 
    min_sat_amount = 5
    sat_amount = 9
    method_type = "Newton-Gauss"
    all_all_errors = []

    X_0 = np.array([0,0,6370,0])

    pos = Positions(pos_guess=X_0, toll=1e-8)

    for sat_amount in range(min_sat_amount, sat_amount+1):
        all_errors = []

        rand_thetas, rand_phis = pos.random_sat_angles(sat_pos_amount, sat_amount)
        for i in range(len(rand_phis)):
            max_error = pos.distance_w_error(rand_thetas[i], rand_phis[i], 1e-8, 1e-8, sat_amount, method_type)
            all_errors.append(max_error)

        print(f'Done calculating the error for satellite amount {sat_amount} with the {method_type} method\n')

        fig, ax = plt.subplots()
        ax.set_xlabel('Satellite amount')
        ax.set_ylabel('Error [km]')
        ax.set_title(f'Error w.r.t {sat_amount} satellites for the {method_type} method')
        bp = ax.boxplot(all_errors, positions=[sat_amount])  
        plt.setp(bp['whiskers'], color='k', linestyle='-')
        plt.setp(bp['fliers'], markersize=3.0)
        fig.savefig(f'figures/sat_num_{sat_amount}.png')

        if sat_amount <= min_sat_amount:
            continue
        all_all_errors.append(all_errors)

    fig, ax = plt.subplots()
    ax.set_xlabel('Satellite amount')
    ax.set_ylabel('Error [km]')
    ax.set_title('Error w.r.t the amount of satellites for the {method_type} method')
    bp = ax.boxplot(all_all_errors, positions=[i for i in range(6, sat_amount+1)])  
    plt.setp(bp['whiskers'], color='k', linestyle='-')
    plt.setp(bp['fliers'], markersize=3.0)
    fig.savefig(f'figures/all_sats.png')

    plt.clf()
    plt.plot([i for i in range(6, sat_amount+1)], [np.max(i) for i in all_all_errors], label="max")
    plt.plot([i for i in range(6, sat_amount+1)], [np.min(i) for i in all_all_errors], label="min")
    plt.plot([i for i in range(6, sat_amount+1)], [np.mean(i) for i in all_all_errors], label="mean")
    plt.xlabel("Satellite amount")
    plt.ylabel("Error [km]")
    plt.title("Different errors w.r.t. satellite amount with Newton-Gauss")
    plt.legend(loc='best')
    plt.xticks(np.arange(6, sat_amount+1, 1))
    plt.savefig('figures/diff_errors')
    plt.show()