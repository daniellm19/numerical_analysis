import matplotlib.pyplot as plt
from numpy import sin, pi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Pendulum:

    G=9.81

    def __init__(self, x: list, n: int, T: float, L: float = 2) -> None:
        """Just the init function

        Args:
            x (list):
                A list with the inital value for the angle and the initial derivative value of that angle
            n (int):
                The number of steps in the Euler's method
            T (float):
                The end of the time interval (0,T)
            L (float, optional): 
                The length of the pendulum rod. Defaults to 2.
        """
        self.x = x
        self.n = n
        self.T = T
        self.L = L
        self.h = self.T/self.n

    def __ydot(self, y: list):
        """The differential system

        Args:
            y (list):
                The value of the angle as well as the derivative value

        Returns:
            tuple: A tuple with the two results from the differntial systems
        """
        return (y[1], -self.G/self.L*sin(y[0]))

    def euler_step(self):
        """The Euler method

        Returns:
            _type_: _description_
        """

        t_list = [round(self.h*i,12) for i in range(n+1) if i!=0]    # List of time values          
        y_list = [] # list of angle values
        x = self.x
        for _ in range(n):
            y=[0,0]
            y[0], y[1] = x[0] + self.h * self.__ydot(x)[0], x[1] + self.h * self.__ydot(x)[1]
            y_list.append(y)
            x = y

        return y_list, t_list

    def plot_euler_step(self):
        y, t = self.euler_step()
        plt.figure(figsize=(9, 3))
        plt.plot(t,y, label=[r"$\theta$", r"$\dot{\theta}$"])
        plt.xlabel("Time [s]")
        plt.ylabel("Y vector [rad]")
        plt.title("Simple pendulum swing w.r.t time")
        plt.legend(loc='best')
        plt.show()

    def simple_simulate(self, dt = 0.1):
        # Source: https://scipython.com/book2/chapter-7-matplotlib/problems/p77/animating-a-pendulum/

        def get_coords(th):
            """Return the (x, y) coordinates of the bob at angle th."""
            return (self.L * np.sin(th), -self.L * np.cos(th))

        y_list, _ = self.euler_step()
        # Initialize the animation plot. Make the aspect ratio equal so it looks right.
        fig = plt.figure()
        ax = fig.add_subplot(aspect='equal')
        # The pendulum rod, in its initial position.
        x0, y0 = get_coords(self.x[0])
        line, = ax.plot([0, x0], [0, y0], lw=3, c='k')
        # The pendulum bob: set zorder so that it is drawn over the pendulum rod.
        bob_radius = 0.08
        circle = ax.add_patch(plt.Circle(get_coords(self.x[0]), bob_radius,
                            fc='r', zorder=3))
        # Set the plot limits so that the pendulum has room to swing!
        ax.set_xlim(-self.L*1.2, self.L*1.2)
        ax.set_ylim(-self.L*1.2, self.L*1.2)

        def animate(i):
            """Update the animation at frame i."""
            x, y= get_coords(y_list[i][0])
            line.set_data([0, x], [0, y])
            circle.set_center((x, y))

        interval = dt * 100
        ani = animation.FuncAnimation(fig, animate, frames=100*len(y_list), repeat=False,
                                    interval=interval)
        plt.show()

if __name__ == "__main__":

    # # Program 1
    # T = 10
    # n = 500
    # y_0 = [pi/2,0]
    # pend = Pendulum(y_0, n, T)
    # pend.plot_euler_step()

    # # Program 2
    # T = 20
    # n = 500
    # y_0 = [pi/12,0]
    # pend = Pendulum(y_0, n, T)
    # pend.plot_euler_step()

    # # Program 3
    # T = 20
    # n = 500
    # y_0 = [pi/12,0]
    # pend = Pendulum(y_0, n, T)
    # pend.simple_simulate()

    # Program 4
    T = 20
    n = 500
    y_0 = [pi/2,0]
    pend = Pendulum(y_0, n, T)
    pend.simple_simulate()