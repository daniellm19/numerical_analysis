from numpy import sin, cos, sqrt, pi
rho = 26570 # kilometers
x = 0
y = 0
z = 6370

def location(theta, phi):
    A = rho * sin(phi) * cos(theta) if abs(rho * sin(phi) * cos(theta)) > 1e-10 else 0
    B = rho * sin(phi) * sin(theta) if abs(rho * sin(phi) * sin(theta)) > 1e-10 else 0
    C = rho * cos(phi) if abs(rho * cos(phi)) > 1e-10 else 0
    distance = sqrt(pow((x - A), 2) + pow((y - B), 2) + pow((z - C), 2))
    t = distance / 299792.458

    print('A =', A, 'km\nB =', B, 'km\nC =', C, 'km\nDistance from north pole =', distance, 'km\nt =', t, 'sec')

location(pi,0.1)