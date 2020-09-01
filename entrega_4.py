# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 14:40:30 2020

@author: Fran
"""

from matplotlib.pylab import *
from scipy.integrate import odeint, solve_ivp
from eulerint import eulerint
import math
import numpy as np

#Datos:
m = 1.0 #Kg
f = 1.0 #Hz
e = 0.2 #chi
w = 2*(math.pi)*f
k = m*(w**2)
c = 2*e*w*m

   
def zp(z, t):
    x = z[0]    
    dx = z[1]
    zdot = np.zeros(2)
    zdot[0] = dx
    zdot[1] = -(c/m)*dx - (k/m)*x
    return zdot

z0 = array([1.0, 1.0]) #x, vx  en condicion inicial

t = linspace(0., 4., 100)

sol = odeint(zp, z0, t)

z_odeint = sol[:,0]

# Resuelve la EDO
z_real = np.exp(-c*t/2)*np.cos(w*t)

sol1 = eulerint(zp, z0, t, Nsubdivisiones=1)
sol2 = eulerint(zp, z0, t, Nsubdivisiones=10)
sol3 = eulerint(zp, z0, t, Nsubdivisiones=100)
z_euler1 = sol1[:,0]
z_euler2 = sol2[:,0]
z_euler3 = sol3[:,0]

title("Aproximaciones oscilador armónico")
plot(t, z_odeint, label="odeint", color="b")
plot(t, z_euler1, "--", label="eulerint Nsubdivisiones = 1", color="g")
plot(t, z_euler2, "--", label="eulerint Nsubdivisiones = 10", color="r")
plot(t, z_euler3, "--", label="eulerint Nsubdivisiones = 100", color="orange")
plot(t, z_real, "--", label="real", color="k", linewidth=2)
legend()
xlabel("Tiempo (t)")
ylabel("X(t)")
savefig("oscilador_armonico.png")
show()


