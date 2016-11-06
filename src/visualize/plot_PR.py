import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import math

R_const=1.0
a=2.0/49.0
b=2.0/21.0
w_pian=0.344

G = -1

Tcr = 0.0729
T = (0.6*Tcr)
		   
alfa_T=(1.0+(0.37464+1.54226*w_pian-0.26992*w_pian*w_pian)*(1.0-math.sqrt(T/Tcr)))*(1.0+(0.37464+1.54226*w_pian-0.26992*w_pian*w_pian)*(1.0-math.sqrt(T/Tcr)))

def func(rho):
    pressure = rho*R_const*T/(1.0-b*rho)-a*alfa_T*rho*rho/(1.0+2.0*b*rho-b*b*rho*rho)
    psi = math.sqrt(2.0*(pressure-rho/3.0)*3/G)
    return psi

print [i for i in range(0,10)]