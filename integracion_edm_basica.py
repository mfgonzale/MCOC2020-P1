# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 19:34:06 2020

@author: Fran
"""
from numpy import *
import numpy as np
from sympy import *
import scipy as sp
from scipy.integrate import odeint
import math

omega = 7.25*(10**-5)   #velocidad angular de la Tierra en rad/s
t = sp.linspace(0, 3600*4, 4*3600 + 1) #np.array([0, 86400, 86400*2]) #24 hrs a segundos y 48 horas a segundos
mt = 5.972*(10**24) #Masa de la Tierra
G = 6.67*(10**-11) #constante gravitatoria universal en (m***3)/kg*(s**2)
r = 700000 + 6371000 #distancia del satelite al centro de la Tierra


def radio(posicion):
    r = (np.dot(posicion,posicion))**0.5
    return r
    
    
    
def satelite(z,t):
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

#Condicion inicial:
vt = 25000/3.6 #m/s
z0 = np.array([r,0,0,0,vt,0])  #x, y, z, vx, vy, vz
 

sol = odeint(satelite,z0,t)


import matplotlib.pyplot as plt

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

plt.figure(1)
plt.plot(x,y)
plt.title("Trayectoria de orbita de un satelite")
plt.xlabel("Recorrido en x (m)")
plt.ylabel("Recorrido en y (m)")
plt.grid(True)

plt.savefig("satelite1.png")
plt.show()

plt.figure(2)
plt.title("Trayectoria en el eje x de la orbita del satelite")
plt.ylabel("Recorrido (m)")
plt.xlabel("Tiempo (s)")
plt.grid(True)
plt.plot(t,x)
plt.savefig("satelite2.png")
plt.show()

plt.figure(3)
plt.title("Trayectoria en el eje y de la orbita del satelite")
plt.ylabel("Recorrido (m)")
plt.xlabel("Tiempo (s)")
plt.grid(True)
plt.plot(t,y)
plt.savefig("satelite3.png")
plt.show()

plt.figure(4)
plt.title("Trayectoria en el eje z de la orbita del satelite")
plt.ylabel("Recorrido (m)")
plt.xlabel("Tiempo (s)")
plt.grid(True)
plt.plot(t,z)
plt.savefig("satelite4.png")
plt.show()

plt.figure(5)
plt.title("Trayectoria del radio r(t) de la orbita del satelite")
plt.ylabel("Recorrido (m)")
plt.xlabel("Tiempo (s)")
plt.grid(True)

plt.plot(t,np.sqrt(x**2 + y**2 + z**2))
#Los siguientes codigos para trazar la superficie de la Tierra y la atmosfera
#sobreescribían la función r(t) por lo que no pude incluirlo en el gráfico
#plt.axvline(x = 6371000.0, color='r')
#plt.axvline(6371000+80000, color='g', linestyle='-', label = "Atmósfera")
plt.savefig("satelite5.png")
plt.show()
  
  
Comentarios:  
    Para conseguir vt se tuvo probar con distintos valores para ver la velocidad minima necesaria para que la órbita del satelite no variara de su forma eliptica, y por ende, no caiga hacia la atmosfera terrestre.  
    Finalmente, vt minimo es 25000 m/s  
