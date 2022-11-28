import matplotlib.pyplot as plt
from numpy import sin, pi

#constants:
g=9.81
L=2

def ydot(y: list):
    return (y[1], -g/L*sin(y[0]))

def euler_step(x: int, n: int, T: int):
    """The Euler method

    Args:
        x (int): _description_
        n (int): _description_
        T (int): _description_

    Returns:
        _type_: _description_
    """
    h = T/n                                                 # The step size
    t_list = [round(h*i,12) for i in range(n+1) if i!=0]    # List of time values          
    y_list = []
    for _ in range(n):
        y=[0,0]
        for j in range(len(y)):
            y[j] = x[j] + h * ydot(x)[j]
        y_list.append(y)
        x = y

    return y_list, t_list

def main():
    T = 10
    n = 500
    y_0 = [pi/2,0]
    y, t = euler_step(y_0, n, T)
    plt.figure(figsize=(9, 3))
    plt.plot(t,y)
    plt.show()

main()
