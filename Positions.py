import numpy as np
from numpy import sin, cos, sqrt, pi, mean, std, linalg as LA
from nptyping import NDArray, Shape, Int
import matplotlib.pyplot as plt
import random
from sympy.utilities.iterables import multiset_permutations

class Positions:
    """
    A class that can calculate using either the Newton or Gauss-Newton method the intersect of 
    a possible receiver given n amount of satellites.

    ...

    Attributes
    ----------
    pos_guess : NDArray[Shape["4"], Int]
        numpy array of size 4 that takes in the x, y, z coordinates and the distance d to the satellite that
        the user guesses where the receiver is located
    toll : float
        the tolerance that the Newton or Gauss-Newton method will use to determine when the root has
        been found

    Methods
    -------
    add_satellite_polar(theta: float, phi: float):
        Takes in the altitude and polar angle of the satellite and stores it in the class
    add_satellite_cart(self, x: float, y: float, z: float, t: float):
        Takes in the x, y, z coordinates of the satellite and the time t it takes to receive a signal from it
    remove_all_sat():
        Removes all of the satellites from the class
    find_intersection(self, method : str = None):
        Finds the intersection given the x0 starting guess, can choose to use either the Newton or Gauss-Newton
        method. Defaults to Gauss-Newton when the amount of satellites is over 4 and defaults to Newton when
        its below or equal to 4.
    random_sat_angles(self, pos_amount: int, sat_amount: int = self.sat_amount):
        Generates {pos_amount} of random satellite positions for {sat_amount} of satellites
    def polar_calc(self, A: float, B: float, C: float):
        Calculates the polar coordinates (altitude and azimuth/polar angle) of the satellite given
        its x, y, and z coordinates
    cartesian_calc(self, phi: float, theta: float) -> dict:
        Calculates the cartesian coordinates of the satellites given the altitude (phi) and azimuth/polar angle (theta)
        angles of the satellite
    get_incorr_angle(self, err: int, angles: list) -> list:
        Creates all of the combinations of a error value being either added or subtracted from
        a list of angle coordinates
    calculate_multiple_cartesian(self, theta: float, phi: float) -> float:
        Gets the cartesian values and the distance from the angles (altitude and azimuth/polar) from
        the satellites
    f(self, x: list, A: list, B: list, C: list, t: list) -> list:
        Takes in a 1x4 vector of x, y, z, and d (t_{recv}-t_{sat}) for the initial value along with the
        cartesian positions of the satellites. 
        Creates a matrix and appends each equation with relevant variable data into the list 
        and returns a {rows}x1 matrix of the equations.
    df(self, x: list, A: list, B: list, C: list, t: list) -> list:
        Calculates each row of jacobi matrix independently (Hard coded). Takes in a 5x1 vector
        of initial conditions and creates a {rows}x{rows} matrix
    newton_mult(self, tol: int, A: list, B: list, C: list, t: list) -> list
        Calculates the intercept that the satellites estimate with the Newton method.
        Takes in inital values along with a tolerance to determine when the root has
        been reached.
    newton_gauss_mult(self, tol: int, A: list, B: list, C: list, t: list) -> list:
        Calculates and returns the estimated intercept of the receivers position using the Newton-Gauss method.
    def distance_w_error(self, thetas: list, phis: list, err: int, tol_err: int, sat_amount: int, method_type: str = "Newton-Gauss") -> float:
        Takes a list of angle coordinates of satellites and adds errors to those coordinates,
        then it calculates the error (using either Newton or Newton-Gauss) and returns 
        the maximum error found
    """

    RHO = 26570                 # distance of satellites from earth [km]
    LIGHT_SPEED = 299792.458    # speed of light [km/s]

    def __init__(self, pos_guess: NDArray[Shape["4"], Int], toll: float) -> None:
        self.x0 = pos_guess
        self.toll = toll
        self.sat_amount = 0
        self.theta = []
        self.phi = []
        self.A = []
        self.B = []
        self.C = []
        self.t = []

    def add_satellite_polar(self, theta: float, phi: float):
        """
        Adds a satellite to the class
        
        Arguments
        ----------
        theta : The altitude angle of the satellite from the +z axis
        phi : The azimuth/polar angle of the satellite from the x- to y-axis
        """
        self.sat_amount += 1
        self.theta.append(theta)
        self.phi.append(phi)
        carts = self.cartesian_calc(theta, phi)
        self.A.append(carts["A"])
        self.B.append(carts["B"])
        self.C.append(carts["C"])
        self.t.append(carts["t"])

    def add_satellite_cart(self, x: float, y: float, z: float, t: float):
        """
        Adds a satellite to the class
        
        Arguments
        ----------
        x : The x-axis coordinate of the satellite
        y : The y-axis coordinate of the satellite
        z : The z-axis coordinate of the satellite
        t : The time that it takes for a signal to travel from the satellite to the receiver
        """
        self.sat_amount += 1
        self.A.append(x)
        self.B.append(y)
        self.C.append(z)
        self.t.append(t)

        phi, theta = self.polar_calc(x,y,z)
        self.phi.append(phi)
        self.theta.append(theta)

    def remove_all_sat(self):
        """
        Removes all of the satellites from the class
        """
        self.sat_amount = 0
        self.theta = []
        self.phi = []
        self.A = []
        self.B = []
        self.C = []
        self.t = []


    def find_intersection(self, method : str = None):
        """
        Finds the intersection using the Newton method if the satellites are less than 5
        but uses the Newton-Gauss method if they're more than 5
        
        Arguments
        ----------
        method (optional) : Can specify the method used (either "Newton" or "Newton-Gauss")
        """
        if method == None:
            if self.sat_amount>4:
                return self.newton_gauss_mult(self.toll, self.A, self.B, self.C, self.t)
            else:
                return self.newton_mult(self.toll, self.A, self.B, self.C, self.t)
        elif method == "Newton":
            if self.sat_amount>4:
                raise NameError("Cannot use the Newton method if the amount of satellites are more than 4")
            return self.newton_mult(self.toll, self.A, self.B, self.C, self.t)
        elif method == "Newton-Gauss":
            return self.newton_gauss_mult(self.toll, self.A, self.B, self.C, self.t)

    def random_sat_angles(self, pos_amount: int, sat_amount: int = None):
        '''
        Generates {pos} random positions for {sat_amount} amount of satellites each
        
        Arguments
        ----------
        pos : The number of random positions to calculate for the satellites
        rows : The amount of satellites that become rows in f(x) and df(x) (Jacobi)
        '''
        if sat_amount == None:
            sat_amount = self.sat_amount

        rand_phis = []
        rand_thetas = []
        for _ in range(0, pos_amount):
            rand_phi = []
            rand_theta = []
            for _ in range(0, sat_amount):
                rand_phi.append(random.uniform(0.0, pi/2))
                rand_theta.append(random.uniform(0.0, 2*pi))
            rand_thetas.append(rand_theta)
            rand_phis.append(rand_phi)
        return rand_thetas, rand_phis

    def polar_calc(self, A: float, B: float, C: float):
        phi = np.arctan(A/B)
        theta = np.arccos(C/np.sqrt(pow(A,2)+pow(B,2)+pow(C,2)))

        self.phi.append(phi)
        self.theta.append(theta)

        return phi, theta

        

    def cartesian_calc(self, phi: float, theta: float) -> dict:
        """
        Calculates the cartesian coordinates of the satellites along with the time it takes
        to reach the satellite from the inital position of the receiver

        Arguments
        ----------
        phi : The altitude of the satellite (between 0 and pi/2 for simplification)
        theta : The polar angle of the satellite (between 0 and 2*pi)
        """
        if (phi>(np.pi/2)) or (0>phi):
            raise ValueError("Phi (the altitude) should be between 0 and pi/2")
        #if (theta>(2*np.pi)) or (0>theta):
        #    raise ValueError("Theta (the polar angle) should be between 0 and 2*pi")

        A = self.RHO * sin(phi) * cos(theta) if abs(self.RHO * sin(phi) * cos(theta)) > 1e-10 else 0
        B = self.RHO * sin(phi) * sin(theta) if abs(self.RHO * sin(phi) * sin(theta)) > 1e-10 else 0
        C = self.RHO * cos(phi) if abs(self.RHO * cos(phi)) > 1e-10 else 0
        distance = sqrt(pow((self.x0[0] - A), 2) + pow((self.x0[1] - B), 2) + pow((self.x0[2] - C), 2))
        t = distance / self.LIGHT_SPEED

        ret_dict = {'A': A, 'B': B, 'C': C, 'distance': distance, 't': t}
        return ret_dict

    def get_incorr_angle(self, err: int, angles: list) -> list:
        """
        Creates all of the combinations of a error value being either added or subtracted from
        a list of angle coordinates
        
        Arguments
        ----------
        err : The error that should be used for the permutation
        angles : The list of angles that should be error-riddled
        """
        all_err_perm = []
        err_list = [err] * len(angles)
        for i in range(0,len(err_list)+1):
            err_copy = err_list.copy()
            for j in range(0,i):
                err_copy[j] = -err_list[j]
            multiset_perm = multiset_permutations(err_copy)
            for i in multiset_perm:
                all_err_perm.append(i)

        # Adds all of the error permutations to the a_list
        all_perms = [[a + b for a, b in zip(all_err_perm[i], angles)] for i in range(len(all_err_perm))]

        return all_perms

    def calculate_multiple_cartesian(self, theta: float, phi: float) -> float: #TODO indicate four float value return types
        """
        Gets the cartesian values and the distance from the angles (altitude and azimuth/polar) from
        the satellites

        Arguments
        ----------
        theta : The altitude angles of the satellites
        phi : The azimuth/polar angles of the satellites
        """
        if (len(theta) != len(phi)) or (len(theta)==0) or (len(phi)==0):
            raise ValueError("The thetas and phis should not be zero and should have the same length")

        A, B, C, t = [], [], [], []
        for i in range(len(phi)):
            values = self.cartesian_calc(phi[i], theta[i])     #Derived from prerceived values
            A.append(values['A'])                        #Vector of distances in plane A[km]
            B.append(values['B'])                        #Vector of distances in plane B[km]
            C.append(values['C'])                        #Vector of distances in plane C[km]
            t.append(values['t'])                        #Vector of time for each sat t[s] derived from correct values
        return A, B, C, t

    def f(self, x: list, A: list, B: list, C: list, t: list) -> list: #TODO have return type be np.array
        '''
        Takes in a 1x4 vector of x, y, z, and d (t_{recv}-t_{sat}) for the initial value along with the
        cartesian positions of the satellites. 
        Creates a matrix and appends each equation with relevant variable data into the list 
        and returns a {rows}x1 matrix of the equations.

        Arguments
        ----------
        x : A 1x4 vector of x, y, z, and d (t_{recv}-t_{sat})
        A : A list of the x-axis coordinates of the satellites
        B : B list of the x-axis coordinates of the satellites
        C : C list of the x-axis coordinates of the satellites
        t : The time it takes to travel from the satellite to the receiver (at light speed)
        rows : The amount of satellites (which translates to equations)
        '''
        funcs = []
        for i in range(self.sat_amount):
            funcs.append(pow((x[0]-A[i]), 2) + pow((x[1]-B[i]),2) + pow((x[2]-C[i]),2) - pow(self.LIGHT_SPEED,2) * pow((t[i]-x[3]), 2))
        return funcs

    def df(self, x: list, A: list, B: list, C: list, t: list) -> list: #TODO have return type be np.array
        '''
        Calculates each row of jacobi matrix independently (Hard coded). Takes in a 5x1 vector
        of initial conditions and creates a {rows}x{rows} matrix
        
        Arguments
        ----------
        x : A 1x4 vector of x, y, z, and d (t_{recv}-t_{sat})
        A : A list of the x-axis coordinates of the satellites
        B : B list of the x-axis coordinates of the satellites
        C : C list of the x-axis coordinates of the satellites
        t : The time it takes to travel from the satellite to the receiver (at light speed)
        rows : The amount of satellites (which translates to equations)
        '''
        eq_list = [0]*self.sat_amount
        for i in range(0, len(eq_list)):
            eq_list[i] = [(2*(x[0]-A[i])), (2*(x[1]-B[i])), (2*(x[2]-C[i])), (2*pow(self.LIGHT_SPEED,2)* (t[i]-x[3]))]
        return np.array(eq_list)

    def newton_mult(self, tol: int, A: list, B: list, C: list, t: list) -> list:
        '''
        Calculates the intercept that the satellites estimate with the Newton method.
        Takes in inital values along with a tolerance to determine when the root has
        been reached.

        Arguments
        ----------
        tol : The tolerance of when the root has been reached
        A : A list of the x-axis coordinates of the satellites
        B : B list of the x-axis coordinates of the satellites
        C : C list of the x-axis coordinates of the satellites
        t : The time it takes to travel from the satellite to the receiver (at light speed)
        '''
        x=self.x0
        oldx=x+2*tol
        while LA.norm(x-oldx, np.inf)>tol:
            oldx=x
            s=-LA.solve(self.df(x, A, B, C, t), self.f(x, A, B, C, t))
            x=x+s
        return(x)

    def newton_gauss_mult(self, tol: int, A: list, B: list, C: list, t: list) -> list:
        """
        Calculates and returns the estimated intercept of the receivers position using the Newton-Gauss method.

        Arguments
        ----------
        tol : The tolerance of when the root has been reached
        A : A list of the x-axis coordinates of the satellites
        B : B list of the x-axis coordinates of the satellites
        C : C list of the x-axis coordinates of the satellites
        t : The time it takes to travel from the satellite to the receiver (at light speed)
        """

        x = self.x0
        oldx=x+2*tol
        while LA.norm(x-oldx, np.inf)>tol:
            oldx = x
            left_side = np.dot(np.transpose(self.df(x, A, B, C, t)), self.df(x, A, B, C, t))
            right_side = np.dot((np.transpose(self.df(x, A, B, C, t))), self.f(x, A, B, C, t))
            s = -LA.solve(left_side, right_side)
            x = x+s
        return(x)
        
    def bisection_error(self, a: float, b: float, theta: list, phi: list):
        """
        Calculates and returns the estimated intercept of the receivers position using the Newton-Gauss method.

        Arguments
        ----------
        a : Lower bound on the error
        b : Upper bound on the error
        theta : The angles 
        """

        if 0.0001 - self.distance_w_error(theta, phi, a) * self.distance_w_error(theta, phi, b) >= 0:
            return None
        else:
            fa=distance_w_error(theta, phi, a)
            while (b-a)/2>self.toll:
                c=(a+b)/2
                fc=self.distance_w_error(theta, phi, c)
                if fc==0:break
                if fc*fa<0:
                    b=c
                else:
                    a=c
                    fa=fc          
        return((a+b)/2)

    def distance_w_error(self, thetas: list, phis: list, incorr_phis: list, method_type: str = "Newton-Gauss") -> float:
        '''
        Takes a list of angle coordinates of satellites and adds errors to those coordinates,
        then it calculates the error (using either Newton or Newton-Gauss) and returns 
        the maximum error found

        Arguments
        ----------
        thetas : The altitude angles
        phis: The azimuth/polar angles
        err : The error you want to add to the satellites
        tol_err : The tolerance error for either the Newton or Newton-Gauss method
        sat_amount : The amount of satellites
        method_type: Specifies either the Newton method ("Newton") or the Newton-Gauss method ("Newton-Gauss")
        '''

        if (method_type!="Newton-Gauss") and (method_type!="Newton"):
            raise ValueError("method_type should be either 'NG' (for Newton-Gauss) or 'N' (for Newton)")

        old_sat_amount = self.sat_amount
        self.sat_amount = sat_amount

        all_errors = []
        for incorr_phi in incorr_phis:
            A, B, C, _ = self.calculate_multiple_cartesian(thetas, incorr_phi)
            _, _, _, t = self.calculate_multiple_cartesian(thetas, phis)
            if method_type=="Newton-Gauss":
                x,y,z,_ = self.newton_gauss_mult(tol_err, A, B, C, t)
            elif method_type=="Newton":
                x,y,z,_ = self.newton_mult(tol_err, A, B, C, t)
            all_errors.append(sqrt(pow(x - self.x0[0], 2) + pow(y - self.x0[1], 2) + pow(z - self.x0[2], 2)))
        
        self.sat_amount = old_sat_amount             # Resets the satellite amount since none are created
        return max(all_errors)

    def distance_w_max_error(self, theta: list, phi: list, incorr_phis: list, allowed_error: float, method_type: str = "Newton-Gauss"):
        return self.distance_w_error(theta, phi, incorr_phis) - allowed_error

    def bisection_error(self, theta: list, phi: list, a: float, b: float, tol: float, allowed_error: float, method_type: str = "Newton-Gauss"):
        '''gert ráð fyrir að búið se að skilgreina f(x) fyrir utan t.d.
        def f(x):
            return(x**2-2)
        '''
        a_incorr_phis = self.get_incorr_phis(a, phi)
        b_incorr_phis = self.get_incorr_phis(b, phi)
        if self.distance_w_max_error(theta, phi, a_incorr_phis, allowed_error)*self.distance_w_max_error(theta, phi, b_incorr_phis, allowed_error) >= 0:
            print("Bisection method failed.")
            return None
        else:
            fa=self.distance_w_max_error(theta, phi, a_incorr_phis, allowed_error)
            while (b-a)/2>tol:
                c=(a+b)/2
                c_incorr_phis = self.get_incorr_phis(c, phi)
                fc=self.distance_w_max_error(theta, phi, c_incorr_phis, allowed_error)
                if fc==0:break
                if fc*fa<0:
                    b=c
                else:
                    a=c
                    fa=fc
        return((a+b)/2)

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