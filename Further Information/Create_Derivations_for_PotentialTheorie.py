"""
Created on Tue Jan 18 10:24:39 2022

@author: ei04yzyc
"""

from sympy import symbols, log, diff, integrate, simplify
from scipy.integrate import dblquad

#Symbols
xi, eta, x, y, z = symbols ("xi eta x y z")

#Abstand
rho = ((xi-x)**2+(eta-y)**2+z**2)**(1/2)

#Hilfsvariable
Omega = z*log(rho+z)-rho

#Ableitungen von Omega
dOmega_dx = diff(Omega, x)
dOmega_dy = diff(Omega, y)
dOmega_dz = diff(Omega, z)

d2Omega_dx2 = diff(dOmega_dx, x)
d2Omega_dy2 = diff(dOmega_dy, y)
d2Omega_dz2 = diff(dOmega_dz, z)
d2Omega_dxdy = diff(dOmega_dx, y)
d2Omega_dxdz = diff(dOmega_dx, z)
d2Omega_dydz = diff(dOmega_dy, z)

d3Omega_dx3 = diff(d2Omega_dx2, x)
d3Omega_dy3 = diff(d2Omega_dy2 ,y)
d3Omega_dz3 = diff(d2Omega_dz2, z)
d3Omega_dx2dy = diff(d2Omega_dx2, y)
d3Omega_dx2dz = diff(d2Omega_dx2, z)
d3Omega_dy2dx = diff(d2Omega_dy2, x)
d3Omega_dy2dz = diff(d2Omega_dy2, z)
d3Omega_dz2dx = diff(d2Omega_dz2, x)
d3Omega_dz2dy = diff(d2Omega_dz2, y)
d3Omega_dxdydz = diff(d2Omega_dxdy, z)

d4Omega_dx4 = diff(d3Omega_dx3, x)
d4Omega_dy4 = diff(d3Omega_dy3, y)
d4Omega_dz4 = diff(d3Omega_dz3, z)
d4Omega_dx3dy = diff(d3Omega_dx3, y)
d4Omega_dx3dz = diff(d3Omega_dx3, z)
d4Omega_dy3dx = diff(d3Omega_dy3, x)
d4Omega_dy3dz = diff(d3Omega_dy3, z)
d4Omega_dz3dx = diff(d3Omega_dz3, x)
d4Omega_dz3dy = diff(d3Omega_dz3, y)
d4Omega_dx2dy2 = diff(d3Omega_dx2dy, y)
d4Omega_dx2dz2 = diff(d3Omega_dx2dz, z)
d4Omega_dy2dz2 = diff(d3Omega_dy2dz, z)
d4Omega_dx2dydz = diff(d3Omega_dx2dy, z)
d4Omega_dy2dxdz = diff(d3Omega_dy2dx, z)
d4Omega_dz2dxdy = diff(d3Omega_dz2dx, y)

#Test auswählen und als Beispiel ausgeben
aktuell = d4Omega_dz3dy

print("======================================")
print("Vollständige Ableitung:")
print(aktuell)
print()
print("Vereinfachte Ableitung:")
print(simplify(aktuell))
print("======================================")
