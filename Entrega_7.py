# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 23:18:41 2020

@author: Fran
"""
from matplotlib.pylab import *
from scipy.integrate import odeint
import datetime as dt
from leer_eof import leer_eof

from sys import argv
#eofname = argv[1]
sat_t, sat_x, sat_y, sat_z, sat_vx, sat_vy, sat_vz = leer_eof("S1B_OPER_AUX_POEORB_OPOD_20200828T111242_V20200807T225942_20200809T005942.EOF")

eofname = ("S1B_OPER_AUX_POEORB_OPOD_20200828T111242_V20200807T225942_20200809T005942.EOF")


correccion = 2

tmax = max(sat_t)

hr = 3600. #s
km = 1e3 #m
km3 = (1e3)**3
km5 = (1e3)**5
km6 = (1e3)**6
mu = 398600.440*km3
J2 = 1.75553e10*km5
J3 = -2.61913e11*km6

Nt = 9000
dt = 10.
Radio = 6371. * km  #km
mt = 5.972e24 #Kg
G = 6.67408e-11 #(m***3)/kg*(s**2)
omega = 7.2921150e-5 #rad/s
H0 = 700. * km

zp = zeros(6)

def zpunto(z, t):
    theta = omega*t
    R = np.array([[cos(theta), -sin(theta), 0], 
    [sin(theta), cos(theta), 0],
    [0, 0, 1]])

    R_prima =  np.array([[-sin(theta), -cos(theta), 0], 
    [cos(theta), -sin(theta), 0],
    [0, 0, 0]]) * omega

    R_prima_prima =  np.array([[-cos(theta), sin(theta), 0], 
    [-sin(theta), -cos(theta), 0],
    [0, 0, 0]]) * (omega**2)
    
    x = z[0:3]
    xp = z[3:6]
    
    r = sqrt(dot(x,x))
    
    x_vector = R@x    #vector 
    rnorm = x_vector / r
    Fg = -(mu)/(r**2) * rnorm
    
    
    z2 = x_vector[2]**2     #termino z**2
    rflat = x_vector[0]**2 + x_vector[1]**2      #r = x**2 + y**2
    FJ2 = J2 * x_vector / (r**7)
    FJ2[0] = FJ2[0] * (6*z2 - 1.5 * rflat)  
    FJ2[1] = FJ2[1] * (6*z2 - 1.5 * rflat)
    FJ2[2] = FJ2[2] * (3*z2 - 4.5 * rflat)
    
    FJ3 = zeros(3)
    FJ3[0] = J3 * x_vector[0]*x_vector[2] / r**9 * (10*z2 - 7.5*rflat)
    FJ3[1] = J3 * x_vector[1]*x_vector[2] / r**9 * (10*z2 - 7.5*rflat)
    FJ3[2] = J3 * 1 / r**9 * (4*z2*(z2 - 3*rflat) + 1.5*(rflat)**2)
    
    
    zp = zeros(6)
    zp[0:3] = xp
    
    if correccion == 0:
        zp[3:6] = R.T@(Fg - (2*R_prima@xp + R_prima_prima@x))
    elif correccion == 1:
        zp[3:6] = R.T@(Fg + FJ2 - (2*R_prima@xp + R_prima_prima@x))
    elif correccion == 2:
        zp[3:6] = R.T@(Fg + FJ2 + FJ3 - (2*R_prima@xp + R_prima_prima@x))
        
        
    return zp

#vector de tiempo
t = sat_t   #linspace(0, tmax, 100001)

#Parte en el origen 
x0 = Radio + H0
vt = 24548/3.6 #m/s

z0 = array([sat_x[0], sat_y[0], sat_z[0], sat_vx[0], sat_vy[0], sat_vz[0]])

sol = odeint(zpunto, z0, t)

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

vx = sol[:,3]
vy = sol[:,4]
vz = sol[:,5]



delta = sqrt((x - sat_x)**2 + (y - sat_y)**2 + (z - sat_z)**2)



figure(1)
subplot(3,1,1)
plot(sat_t/3600.,sat_x/1000)
plot(t/3600.,x/1000)
ylabel("$X$ (KM)")
subplot(3,1,2)
plot(sat_t/3600.,sat_y/1000)
plot(t/3600.,y/1000)
ylabel("$Y$ (KM)")
subplot(3,1,3)
plot(sat_t/3600.,sat_z/1000, label="real")
plot(t/3600.,z/1000, label="predicha")
ylabel("$Z$ (KM)")
suptitle("Posicion")
xlabel("Tiempo, $t$ (horas)")
tight_layout(rect=[0,0.03,1,0.95])
savefig("posiciones.png")
      

show()
figure(2) 
plot(t/3600,delta/1000)
suptitle(f"Distancia entre posicion real y predicha")
ylabel("Deriva, $\\delta$ (KM)")
xlabel("Tiempo, $t$ (horas)")
tight_layout(rect=[0,0.03,1,0.95])
savefig("Deriva.png")

show()



import xml
import xml.etree.ElementTree as ET

eofname = ("S1B_OPER_AUX_POEORB_OPOD_20200828T111242_V20200807T225942_20200809T005942.EOF")


def eof_a_pred(eofname, x, y, z, vx, vy, vz):
    
   
    tree = ET.parse(eofname)
    root = tree.getroot()
    
    Data_Block = root.find("Data_Block")		
    List_of_OSVs = Data_Block.find("List_of_OSVs")
    
    for i, osv in enumerate(List_of_OSVs):
    		UTC = osv.find("UTC").text[4:]
    		
    		osv.find("X").text = str(x[i])  #conversion de string a double es implicita
    		osv.find("Y").text = str(y[i])
    		osv.find("Z").text = str(z[i])
    		osv.find("VX").text = str(vx[i])
    		osv.find("VY").text = str(vy[i])
    		osv.find("VZ").text = str(vz[i])
            
    tree.write('S1B_OPER_AUX_POEORB_OPOD_20200828T111242_V20200807T225942_20200809T005942.PRED')
    
    return tree
        
print(eof_a_pred(eofname, x, y, z, vx, vy, vz))



