import math

def f(x):
    return pow(x,3)-4*pow(x,2)+3*x+2

def Df(x):
    return 3*pow(x,2)-8*x+3

def newton(x0,tol):
    oldx=x0+2*tol
    x=x0
    counter = 0
    while abs(oldx-x)>tol:
        oldx=x
        counter += 1
        x=x-f(x)/Df(x)
    return(x)

def newton_searching(tol, leap, start_point, end_point):
    pos = start_point
    root_list = []
    while(pos<end_point):
        print(pos)
        num = newton(pos,tol)
        counter = 0
        for i in root_list:
            if i==num:
                counter = 1
        if counter != 1:
            root_list.append(num)
        pos += leap
    return root_list
        

if "__main__" == __name__:
    newton(-100,0.0000001)
    root_list = newton_searching(0.0000001, 0.5, -200, 200)
    for i in root_list:
        print(i)