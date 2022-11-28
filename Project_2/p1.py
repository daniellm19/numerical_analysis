import matplotlib.pyplot as plt
from numpy import sin, pi
import matplotlib.animation as animation

#constants:
g=9.81
L=2



def ydot(t:float, y: list):
    z = []
    z.append(y[1])
    z.append(-g/L*sin(y[0]))
    return z

def eulerstep(x, n, T):
    h = T/n
    t = 0
    t_list = []
    y_list = []
    for i in range(n):
        y=[0,0]
        for j in range(len(y)):
            y[j]=x[j] + h * ydot(t,x)[j]
        y_list.append(y)
        x = y
        t += h
        t_list.append(t)

    return y_list, t_list

def main():
    T = 10
    n = 500
    y_0 = [pi/2,0]
    y, t = eulerstep(y_0, n, T)
    plt.figure(figsize=(9, 3))
    plt.subplot(121)
    plt.show()




main()
