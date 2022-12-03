import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos, pi
import matplotlib.animation as animation
from collections import deque

#Global constants
g = 9.81

def error(k):
    return pow(10,-k)

def ydot(t:float, y: list, L: float, m:float):
    theta_1pp = ((m*L*pow(y[1],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m*g*sin(y[2])*cos(y[2] - y[0])) + m * L * pow(y[3],2) * sin(y[2] - y[0]) - (m+m) * g*sin(y[0]))/((m+m)*L - m*L*pow(cos(y[2] - y[0]),2))
    theta_2pp = (-1*(m*L*pow(y[3],2)*sin(y[2] - y[0])*cos(y[2] - y[0])) + (m+m)*(g*sin(y[0])*cos(y[2] - y[0])-L*pow(y[1], 2) * sin(y[2] - y[0])-g*sin(y[2])))/((m+m)*L - m*L*pow(cos(y[2] - y[0]), 2))
    return np.array([y[1], theta_1pp, y[3], theta_2pp])


def runge_kutta(x, n, T, L: float, m):
    h = T/n
    t = 0
    t_list = []
    y_list = []
    for _ in range(n):
        k1 = ydot(t, x, L, m)
        k2 = ydot(t + h/2, x + h/2 * k1, L, m)
        k3 = ydot(t + h/2, x + h/2 * k2, L, m)
        k4 = ydot(t + h, x + h * k3, L, m)
        y = np.array(x + h * (k1/6 + k2/3 + k3/3 + k4/6))
        y_list.append(y)
        x = y
        t += h
        t_list.append(t)
    return y_list, t_list, h

def plot(angle: list, angle_error: list, t, title):
    fig = plt.figure(figsize=(12, 7))  
    plt.subplots_adjust(top=0.9, left=0.06, right=0.96, hspace=0.33, bottom=0.08)
    
    for i in range(len(angle_error)):
        plt.subplot(3, 2, i+1)
        plt.plot(t, angle_error[i], label=f"Pendulum's angle (with error 10^-{i+1}) [rad]")
        plt.plot(t, angle, label="Pendulum's angle (no error) [rad]")
        
        plt.xlabel('Time [s]')
        plt.ylabel('Radians')
        plt.legend(loc=2, prop={'size': 6})
        fig.suptitle(title, fontsize=15)
        plt.grid()
        
def error_list(yn: list, y: list, n):
    ret_arr = []
    for i in range(n):
        ret_arr.append([abs(y[i][0]-yn[i][0]),abs(y[i][1]-yn[i][1]),abs(y[i][2]-yn[i][2]),abs(y[i][3]-yn[i][3])])
    return ret_arr

def print_the_error(error: list, t):
    fig = plt.figure(figsize=(12, 7))  
    plt.subplots_adjust(top=0.9, left=0.06, right=0.96, hspace=0.5, bottom=0.08)
    fig.suptitle("Difference in angle with starting error 10^(-k) with k ∈ {1,2,3,4,5}")
    counter = 1
    for element in error:
        plt.subplot(3, 2, counter)
        plt.grid()
        angle1, velocity1, angle2, velocity2 = map(list, zip(*element))
        plt.plot(t, angle1, label="Inner pendulum error")
        plt.plot(t, angle2, label="Outer pendulum error")
        plt.xlabel('Time [s]',fontsize=7)
        plt.ylabel('Radians', fontsize=7)
        plt.legend(loc=2, prop={'size': 6})
        plt.title(f"Starting with error 10^-({counter}) [rad]", fontsize=9)
        counter += 1 
    plt.show()
        
def animate_penduli(x_1, y_1, x_2, y_2, ex_1, ey_1, ex_2, ey_2, n, h, k):
    
    fig = plt.figure(figsize=(10, 5))
    plt.subplots_adjust(top=0.9, left=0.06, right=0.96, hspace=0.3, bottom=0.06)
    ax = fig.add_subplot(1, 2, 1, autoscale_on=False, xlim=(-4.2, 4.2), ylim=(-4.2, 4.2))
    ax.set_title("Double pendulum")
    ax2 = fig.add_subplot(1, 2, 2, autoscale_on=False, xlim=(-4.2, 4.2), ylim=(-4.2, 4.2))
    ax2.set_title(f"Double pendulum with error 10^-{k}")
    ax.grid()   
    ax2.grid()
    
    # inital setup for double pendulum
    line_1, = ax.plot([], [], 'o-', c='blue', lw=1.5)
    line_2, = ax.plot([], [], 'o-', c='red', lw=1.5)
    trace, = ax.plot([], [], '.-', c='red', lw=0.5, ms=1)
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    trajectory_x, trajectory_y = deque(maxlen=n), deque(maxlen=n)
    
    # inital setup for double pendulum with error k
    eline_1, = ax2.plot([], [], 'o-', c='blue', lw=1.5)
    eline_2, = ax2.plot([], [], 'o-', c='red', lw=1.5)
    etrace, = ax2.plot([], [], '.-', c='red', lw=0.5, ms=1)
    etime_text = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)
    etrajectory_x, etrajectory_y = deque(maxlen=n), deque(maxlen=n)
    

    
    def animate(i):
        # animating double pendulum
        x1 = [0, x_1[i]]
        y1 = [0, y_1[i]]
        x2 = [x_1[i], x_2[i]]
        y2 = [y_1[i], y_2[i]]
        
        # animating double pendulum with error k
        ex1 = [0, ex_1[i]]
        ey1 = [0, ey_1[i]]
        ex2 = [ex_1[i], ex_2[i]]
        ey2 = [ey_1[i], ey_2[i]]
        
        if i == 0:
            trajectory_x.clear()
            trajectory_y.clear()
            etrajectory_x.clear()
            etrajectory_y.clear()

        trajectory_x.appendleft(x2[1])
        trajectory_y.appendleft(y2[1])
        etrajectory_x.appendleft(ex2[1])
        etrajectory_y.appendleft(ey2[1])

        
        line_1.set_data(x1, y1)
        line_2.set_data(x2, y2)
        trace.set_data(trajectory_x, trajectory_y)
        time_text.set_text(f"time = {i*h:.1f}s")
        
        eline_1.set_data(ex1, ey1)
        eline_2.set_data(ex2, ey2)
        etrace.set_data(etrajectory_x, etrajectory_y)
        etime_text.set_text(f"time = {i*h:.1f}s")
        
        return line_1, line_2, eline_1, eline_2, trace, etrace, time_text, etime_text
   
    ani = animation.FuncAnimation(
        fig, animate, len(x_1), interval=h*1000, blit=True, repeat=False)
    plt.show()
    return 

# error table for difference in angle at t = 40 (last value), for different k
def make_table(error: list):
    print("Error at T=40")
    print("\tInner\t\t      Outer")
    counter = 1
    for item in error:
        print(f'K = {counter}:  {item[-1][0]}   {item[-1][2]}')
        counter += 1
    print()
def main():
    T = 40
    n = 500
    L = 2
    m = 1
    k = [1,2,3,4,5]
    
    # calculating runge kutta for double pendulum with no error
    y0 = np.array([2*pi/3, 0, pi/6, 0])
    y, t, h = runge_kutta(y0, n, T, L, m)
    angle1, velocity1, angle2, velocity2 = map(list, zip(*y))
    
    # calculating runge kutta and errors for all values of k
    diff_y = []
    angle1_error_list = []
    angle2_error_list = []
    error_difference = []
    
    
    for num in k:
        y0 = np.array([2*pi/3+error(num), 0, pi/6 + error(num), 0])
        y_n, t_n, h_n = runge_kutta(y0, n, T, L, m)
        diff_y.append([y_n, t_n, h_n])
        angle1_error, velocity1_error, angle2_error, velocity2_error = map(list, zip(*y_n))
        angle1_error_list.append(angle1_error), angle2_error_list.append(angle2_error)
        y_n_error = error_list(y_n, y, n)
        error_difference.append(y_n_error)
        
    # position comparison graph
    if input("Do you want to see position comparison, shown on a graph? (y/n) ").lower() == 'y':
        plot(angle1, angle1_error_list, t, 'Inner pendulum comparison')
        plt.show()
        plot(angle2, angle2_error_list, t, 'Outer pendulum comparison') 
        plt.show()
        
     # error comparison graph
    if input("Do you want to see the difference in angle, with starting error 10^(-k) where k = {1, 2, 3, 4, 5}, shown on a graph, as well as getting a table of error difference at T = 40? (y/n) ") == 'y':
        print_the_error(error_difference, t)
        make_table(error_difference)
        
    # position comparison animation
    if input("Do you want to see the position comparison, with chosen k value, animated? (y/n) ") == 'y':
        compare_k = input('Animate correct double penduli, and compare with the one with starting error 10^(-k) where k = {1, 2, 3, 4, 5}. Choose k: ')
        if not compare_k.isdigit():
            print("wrong input, I´m shutting down!")
            quit()
        k = int(compare_k)
        if k < 1 | k > 5:
            print("wrong input, I´m shutting down!")
            quit()
        
        x_1, y_1 = L * sin(angle1[:]), -L * cos(angle1[:])
        x_2, y_2 = L * sin(angle2[:]) + x_1, -L * cos(angle2[:]) + y_1
        ex_1, ey_1 = L * sin(angle1_error_list[k-1][:]), -L * cos(angle1_error_list[k-1][:])
        ex_2, ey_2 = L * sin(angle2_error_list[k-1][:]) + ex_1, -L * cos(angle2_error_list[k-1][:]) + ey_1
        animate_penduli(x_1, y_1, x_2, y_2, ex_1, ey_1, ex_2, ey_2, n, h, k)
    

if __name__ == '__main__':
    main()