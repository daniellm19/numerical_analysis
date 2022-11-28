import matplotlib.pyplot as mp

def f(t, y):
    return t*y-pow(t, 3)

def code(y0, n, T):
    h = T/n
    y = [y0]
    t = [0] 
    for i in range(n):
        t.append(t[i]+h)
        y.append(y[-1] + f(t[-1], y[-1]) * h)
    mp.plot(t, y)

def answer(t):
    return pow(t, 2) + 2

def main():
    code(2, 50, 3)
    h = 0.01
    y = [answer(0)]
    t = [0]
    while t[-1] <= 3:
        y.append(answer(t[-1]))
        t.append(t[-1] + h)
    mp.plot(t, y)
    mp.show()
    
    return 
main()
