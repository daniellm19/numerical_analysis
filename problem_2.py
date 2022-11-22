from numpy import sin, cos, sqrt, pi
rho = 26570 # kilometers
x = 0
y = 0
z = 6370

def location(phi, theta):
    A = rho * sin(phi) * cos(theta) if abs(rho * sin(phi) * cos(theta)) > 1e-10 else 0
    B = rho * sin(phi) * sin(theta) if abs(rho * sin(phi) * sin(theta)) > 1e-10 else 0
    C = rho * cos(phi) if abs(rho * cos(phi)) > 1e-10 else 0
    distance = sqrt(pow((x - A), 2) + pow((y - B), 2) + pow((z - C), 2))
    t = distance / 299792.458

    print('A =', A, 'km\nB =', B, 'km\nC =', C, 'km\nDistance from north pole =', distance, 'km\nt =', t, 'sec\n\n')
    
    return

#Problem 3 starts here:
corr_phi = [pi/8, pi/6, 3*pi/8, pi/4]
corr_theta = [-pi/4, pi/2, 2*pi/3, pi/6]

def main():
    for i in range(len(corr_phi)):
        location(corr_phi[i], corr_theta[i])
    return

if __name__ == "__main__":
    main()