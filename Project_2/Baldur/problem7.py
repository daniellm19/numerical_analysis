import matplotlib as plt
import numpy as np
from numpy import sin,cos

#Global constants
g = 9.81

def theta1_twoprime(L1,L2,m1,m2,theta1_prime,theta2_prime,theta1,theta2):
    delta = theta2-theta1
    numerator = (m2*L1*pow(theta1_prime,2)*sin(delta)*cos(delta)) + (m2*g*sin(theta2)*cos(delta)) + m2 * L2 * pow(theta2_prime,2) * sin(delta)-((m1+m2)g*sin(theta1)) 
    denominator = (m1+m2)*L1 - m2*L1*pow(cos(delta),2)
    return numerator/denominator

def theta2_twoprime(L1,L2,m1,m2,theta1_prime,theta2_prime,theta1,theta2):
    delta = theta2-theta1
    numerator = -1*(m2*L1*pow(theta2_prime,2)*sin(delta)*cos(delta)) + (m1+m2)*(g*sin(theta1)*cos(delta)-L1*pow(theta1_prime, 2)*sin(delta)-g*sin(theta2))
    denominator = (m1+m2)*L2 - m2*L2*pow(cos(delta), 2)
    return numerator/denominator

