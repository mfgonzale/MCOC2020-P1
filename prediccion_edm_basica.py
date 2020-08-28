# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 13:52:51 2020

@author: Fran
"""

from matplotlib.pylab import *
from scipy.integrate import odeint
import numpy as np
import math
from sympy import *
import scipy as sp
from numpy import *

hr = 3600. #s
km = 1e3 #m
Radio = 6371. * km  #km
mt = 5.972e24 #Kg
G = 6.67408e-11 #(m***3)/kg*(s**2)
omega = 7.2921150e-5 #rad/s
H0 = 700. * km

#Fgmax = G*Mt/Radio**2

def radio(posicion):
    r = (np.dot(posicion,posicion))**0.5
    return r
    
    
    
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
    
    a = -np.transpose(R)
    zp = np.zeros(6)
    zp[0:3] = z[3:6]
    zp[3:6] = np.dot((-G*mt/((radio(z[0:3]))**3)),z[0:3]) + np.dot(a,((np.dot(R_prima_prima,z[0:3])) + 2*np.dot(R_prima,z[3:6])))  #z2 derivado
    return zp


from datetime import datetime
ti = "2020-08-01T23:00:19.000000"
ti = ti.split("T")
ti = "{} {}".format(ti[0],ti[1])
ti = datetime.strptime(ti, '%Y-%m-%d %H:%M:%S.%f')

tf = "2020-08-03T01:00:19.000000"
tf = tf.split("T")
tf = "{} {}".format(tf[0],tf[1])
tf = datetime.strptime(tf, '%Y-%m-%d %H:%M:%S.%f')

deltaT = (tf-ti).total_seconds()



#Inicial
"""<TAI>TAI=2020-08-01T23:00:19.000000</TAI>
<UTC>UTC=2020-08-01T22:59:42.000000</UTC>
<UT1>UT1=2020-08-01T22:59:41.792858</UT1>
<Absolute_Orbit>+33719</Absolute_Orbit>
<X unit="m">1302228.553340</X>
<Y unit="m">-998548.249232</Y>
<Z unit="m">-6889468.798306</Z>
<VX unit="m/s">1688.511983</VX>
<VY unit="m/s">-7246.635365</VY>
<VZ unit="m/s">1369.824542</VZ> """

#Final
"""<TAI>TAI=2020-08-03T01:00:19.000000</TAI>
<UTC>UTC=2020-08-03T00:59:42.000000</UTC>
<UT1>UT1=2020-08-03T00:59:41.793917</UT1>
<Absolute_Orbit>+33735</Absolute_Orbit>
<X unit="m">2037190.654132</X>
<Y unit="m">5911917.443308</Y>
<Z unit="m">-3321169.495775</Z>
<VX unit="m/s">443.390652</VX>
<VY unit="m/s">-3820.775138</VY>
<VZ unit="m/s">-6540.429894</VZ>"""


x_i = 1302228.553340
y_i = -998548.249232
z_i = -6889468.798306

vx_i = 1688.511983
vy_i = -7246.635365
vz_i = 1369.824542


x_f = 2037190.654132
y_f = 5911917.443308
z_f = -3321169.495775

vx_f = 443.390652
vy_f = -3820.775138
vz_f = -6540.429894

#Vector de tiempo

# t = linspace(0, 5.*hr, 1001)
t = np.linspace(0, deltaT, 9361)

#Parte en el origen
x0 = Radio + H0
vt = 24548/3.6 #m/s
z0 = np.array([x_i, y_i, z_i, vx_i, vy_i, vz_i])

sol = odeint(zpunto, z0, t)


x = sol[:,0:3]

pos_final = np.array([x_f, y_f, z_f, vx_f, vy_f, vz_f]) - sol[-1]

for el in pos_final:
    print(el)
    
"""H = sqrt(x[:, 0]**2 + x[:, 1]**2 + x[:, 2]**2) - Radio

figure()

for i in range(3):
    subplot(3, 1, 1 + i)
    grid(True)
    plot(t / hr, x[:, i])
    
figure()
grid(True)
plot(t / hr, H / km)
axhline(80., linestyle="- -")
axhline(0., linestyle="- -", linewidth=2)

figure()
grid(True)
plot(x[:, 0], x[:, 1])

th = linspace(0, 2 * pi, 400)

axis("equal")

fig = plt.figure()
ax.plt.axes(projection="3d")
ax.plot3D(x[:, 0], x[:, 1], x[:, 2])
show()"""