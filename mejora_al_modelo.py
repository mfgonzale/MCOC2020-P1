# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 18:02:16 2020

@author: Fran
"""

from matplotlib.pylab import *
from leer_eof import leer_eof


t, x, y, z, vx, vy, vz = leer_eof("S1B_OPER_AUX_POEORB_OPOD_20200828T111242_V20200807T225942_20200809T005942.EOF")


figure()

subplot(3, 1, 1)
plot(t,x)
ylabel("X (KM)")
subplot(3, 1, 2)
plot(t,y)
ylabel("Y (KM)")
subplot(3, 1, 3)
plot(t,z)
ylabel("Z (KM)")
xlabel("Tiempo, t (horas)")

show()



#Vector de condicion inicial:
    
z0 = array([x[0], y[0], z[0], vx[0], vy[0], vz[0]])
zf = array([x[-1], y[-1], z[-1], vx[-1], vy[-1], vz[-1]])

#Datos:
hr = 3600. #s
km = 1e3 #m
r = (6371 + 700)*1000 #m
Radio = 6371. * km  #km
mt = 5.972e24 #Kg
G = 6.67408e-11 #(m***3)/kg*(s**2)
omega = 7.2921150e-5 #rad/s
H0 = 700. * km

#Fgmax = G*Mt/Radio**2
import math
import numpy as np
import scipy as sp
from scipy.integrate import odeint
from eulerint import eulerint        
 
def zpunto(z,t):
    theta = omega*t
    R = np.array([[math.cos(theta), -math.sin(theta), 0], 
    [math.sin(theta), math.cos(theta), 0],
    [0, 0, 1]])

    R_prima =  np.array([[-math.sin(theta), -math.cos(theta), 0], 
    [math.cos(theta), -math.sin(theta), 0],
    [0, 0, 0]]) * omega

    R_prima_prima =  np.array([[-math.cos(theta), math.sin(theta), 0], 
    [-math.sin(theta), -math.cos(theta), 0],
    [0, 0, 0]]) * (omega**2)
    
    z1 = z[0:3]

    zp = sp.zeros(6)
    zp[0:3] = z[3:6]
    z1 = (-G*mt/r**3)*z[0:3]  - R.T@(R_prima_prima@z[0:3] + 2*R_prima@z[3:6])
    zp[3:6] = z1
    return zp


sol = odeint(zpunto, z0, t)
z_odeint = sol[:,0]
z_odeint2 = sol[:,1]
z_odeint3 = sol[:,2]
sol1 = eulerint(zpunto, z0, t, Nsubdivisiones=1)
z_euler1 = sol1[:,0]
z_euler2 = sol1[:,1]
z_euler3 = sol1[:,2]

figure()
title("Distancia entre posicion real y predicha")
subplot(2,1,1)
suptitle("Odeint")
plot(t/3600, (z_odeint - z0[0])/1000, label="odeint")
ylabel("Deriva, d (KM)")
subplot(2,1,2)
suptitle("Eulerint")
plot(t/3600, (z_euler1 - z[0])/1000, label="eulerint Nsubdivisiones = 1")
ylabel("Deriva, d (KM)")
xlabel("Tiempo, t (horas)")
show()


