# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:33:23 2020

@author: Christoph Bienefeld
"""

import numpy as np
from scipy.integrate import simps #Alternativ: trapz
from scipy import optimize
import time 


def vergleichsspannung(koords, E, nu, q_z, q_x, q_y, xi_linsp, eta_linsp, xi_vec, eta_vec):
    """Berechnet die Vergleichsspannung im Materialinneren an der mit der Variable koords übergebenen Position"""
    "Die Berechnung der Vergleichsspannung im Materialinneren ist orientiert an: Nikas, Boussinesq-Cerruti functions and a simple technique"
    "Alle benötigten Ableitungen von Omega sind eigenständig mithilfe von Python und der Bibliothek SymPy hergeleitet worden, da diese bei Nikas teilweise fehlerhaft sind"
    
    #Extrahieren der einzelnen Koordinaten
    x, y, z = koords
    
    #Randbedingungen: Formeln für die spätere Integration:
    #Abstand
    # rho = lambda xi,eta: ((xi-x)**2+(eta-y)**2+z**2)**(1/2)
    #Hilfsvariable
    # Omega = lambda xi,eta: z*np.log(rho(xi,eta)+z)-rho(xi,eta)
    #Sowohl rho als auch Omega werden für die Berechnung nicht in Ihrer ursprünglichen Form benötigt
    
    
    #Relevante Ableitungen von Omega
    d3Omega_dx3 = lambda xi,eta: (x - xi)*(z**2 + (eta - y)**2 + (x - xi)**2)**(-12.5)*(3.0*z*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0 + 3.0*z*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 + 2.0*z*(x - xi)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 - 3.0*z*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 - 3.0*z*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 - 3.0*(-x + xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0 + 3.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3
    d3Omega_dy3 = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-12.5)*(-3.0*z*(eta - y)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0 - 3.0*z*(eta - y)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 - 2.0*z*(eta - y)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 3.0*z*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 3.0*z*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 + 3.0*(-eta + y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 3.0*(eta - y)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3
    d3Omega_dz3 = lambda xi,eta: -1.0*z*(1.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**13.5 + 3.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**14.0 + 3.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**14.5 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.0)/(1.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.0 + 3.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.5 + 3.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.0 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.5)
    d3Omega_dx2dy = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-12.5)*(-3.0*z*(eta - y)*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0 - 3.0*z*(eta - y)*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 - 2.0*z*(eta - y)*(x - xi)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 1.0*z*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 1.0*z*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 + 1.0*(-eta + y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 3.0*(eta - y)*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3
    d3Omega_dx2dz = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-13.0)*(3.0*z**2*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 - 3.0*z*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 + z*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(3.0*z + 1.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 1.0*z*(x - xi)**2*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 + 1.0*z*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 - 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (x - xi)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 + 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**12.5 - 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + (x - xi)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**12.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3
    d3Omega_dy2dx = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-12.5)*(3.0*z*(eta - y)**2*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0 + 3.0*z*(eta - y)**2*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 + 2.0*z*(eta - y)**2*(x - xi)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 - 1.0*z*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 - 1.0*z*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 + 3.0*(eta - y)**2*(-x + xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0 + 1.0*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3
    d3Omega_dy2dz = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-13.0)*(3.0*z**2*(eta - y)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 - 3.0*z*(eta - y)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**10.5 + z*(eta - y)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(3.0*z + 1.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.0 + 1.0*z*(eta - y)**2*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 + 1.0*z*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 - 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**11.5 + 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**12.5 - 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + (eta - y)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**12.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3
    d3Omega_dz2dx = lambda xi,eta: (-1.0*x*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.0 - 3.0*x*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.5 - 3.0*x*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.0 - 1.0*x*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.5 + 1.0*xi*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.0 + 3.0*xi*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.5 + 3.0*xi*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.0 + 1.0*xi*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.5)/(1.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.5 + 3.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**17.0 + 3.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**17.5 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**18.0)
    d3Omega_dz2dy = lambda xi,eta: (1.0*eta*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.0 + 3.0*eta*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.5 + 3.0*eta*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.0 + 1.0*eta*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.5 - 1.0*y*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.0 - 3.0*y*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**15.5 - 3.0*y*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.0 - 1.0*y*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.5)/(1.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**16.5 + 3.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**17.0 + 3.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**17.5 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**18.0)
    d3Omega_dxdydz = lambda xi,eta: (eta - y)*(x - xi)*(z**2 + (eta - y)**2 + (x - xi)**2)**(-11.0)*(-3.0*z**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**8.5 + 3.0*z*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**8.5 - z*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(3.0*z + 1.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**9.0 - 1.0*z*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**9.5 + 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**9.5 + 1.0*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**10.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3
    
    d4Omega_dz4 = lambda xi,eta: 1.0*(3.0*z**6*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.0 + 12.0*z**5*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 + 17.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 + 8.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 3.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 - 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 - 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0)/(1.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 + 4.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0 + 6.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.5 + 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.0 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.5)
    d4Omega_dx3dz = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-25.5)*(15.0*z**2*(-x + xi)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 + 15.0*z*(x - xi)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 - z*(x - xi)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(12.0*z + 6.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 - 2.0*z*(x - xi)**3*(3.0*z + 3*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5 + z*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**22.5*(8.0*z*(-x + xi) - 4.0*z*(x - xi) - 3.0*(x - xi)*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)) - 9.0*z*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + 3.0*(-x + xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0 + 3.0*(-x + xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**24.5 + (x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(9.0*z**2 + 3.0*(x - xi)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + (x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5*(6.0*z**2 + 3.0*z*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 3.0*(x - xi)**2) + (x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(3.0*z*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 2.0*(x - xi)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4
    d4Omega_dy3dz = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-25.5)*(15.0*z**2*(eta - y)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 - 15.0*z*(eta - y)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 + z*(eta - y)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(15.0*z + 3.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**22.5 + z*(eta - y)**3*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(12.0*z + 6.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + 2.0*z*(eta - y)**3*(3.0*z + 3*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5 + 9.0*z*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + (-eta + y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(9.0*z**2 + 3.0*(eta - y)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + 3.0*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0 + 3.0*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**24.5 - (eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(3.0*z*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 2.0*(eta - y)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0 + (z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5*(6.0*z**2*(-eta + y) - 3.0*z*(eta - y)*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 2.0*(-eta + y)**3 - 1.0*(eta - y)**3))/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4
    d4Omega_dz3dx = lambda xi,eta: 1.0*z*(3.0*x*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 + 12.0*x*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 + 18.0*x*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 + 12.0*x*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 + 3.0*x*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0 - 3.0*xi*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 - 12.0*xi*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 18.0*xi*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 - 12.0*xi*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 - 3.0*xi*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0)/(1.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.5 + 4.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.0 + 6.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.5 + 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**39.0 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**39.5)
    d4Omega_dz3dy = lambda xi,eta: 1.0*z*(-3.0*eta*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 - 12.0*eta*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 18.0*eta*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 - 12.0*eta*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 - 3.0*eta*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0 + 3.0*y*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 + 12.0*y*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 + 18.0*y*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 + 12.0*y*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 + 3.0*y*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0)/(1.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.5 + 4.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.0 + 6.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.5 + 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**39.0 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**39.5)
    d4Omega_dx2dz2 = lambda xi,eta: 1.0*(3.0*x**2*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**33.5 + 12.0*x**2*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.0 + 18.0*x**2*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 + 12.0*x**2*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 + 3.0*x**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 6.0*x*xi*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**33.5 - 24.0*x*xi*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.0 - 36.0*x*xi*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 - 24.0*x*xi*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 - 6.0*x*xi*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 + 3.0*xi**2*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**33.5 + 12.0*xi**2*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.0 + 18.0*xi**2*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 + 12.0*xi**2*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 + 3.0*xi**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 1.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 - 4.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 - 6.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 - 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5)/(1.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 + 4.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 + 6.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0 + 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.5 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.0)
    d4Omega_dy2dz2 = lambda xi,eta: 1.0*(3.0*eta**2*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**33.5 + 12.0*eta**2*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.0 + 18.0*eta**2*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 + 12.0*eta**2*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 + 3.0*eta**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 6.0*eta*y*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**33.5 - 24.0*eta*y*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.0 - 36.0*eta*y*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 - 24.0*eta*y*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 - 6.0*eta*y*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 + 3.0*y**2*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**33.5 + 12.0*y**2*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.0 + 18.0*y**2*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 + 12.0*y**2*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 + 3.0*y**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 1.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**34.5 - 4.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.0 - 6.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**35.5 - 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 - 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5)/(1.0*z**4*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.0 + 4.0*z**3*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**36.5 + 6.0*z**2*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.0 + 4.0*z*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**37.5 + 1.0*(eta**2 - 2*eta*y + x**2 - 2*x*xi + xi**2 + y**2 + z**2)**38.0)
    d4Omega_dx2dydz = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-25.5)*(15.0*z**2*(eta - y)*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 - 15.0*z*(eta - y)*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 + z*(eta - y)*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(15.0*z + 3.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**22.5 + z*(eta - y)*(x - xi)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(12.0*z + 6.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + 2.0*z*(eta - y)*(x - xi)**2*(3.0*z + 3*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5 + 3.0*z*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + 3.0*(-eta + y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (x - xi)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + 1.0*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0 + 1.0*(eta - y)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**24.5 + (z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5*(2.0*z**2*(-eta + y) - 1.0*z*(eta - y)*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 3.0*(-eta + y)*(x - xi)**2) + (z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(-1.0*z*(eta - y)*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 2.0*(-eta + y)*(x - xi)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4
    d4Omega_dy2dxdz = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-25.5)*(15.0*z**2*(eta - y)**2*(-x + xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 + 15.0*z*(eta - y)**2*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**22.0 - 2.0*z*(eta - y)**2*(x - xi)*(3.0*z + 3*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5 + z*(eta - y)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(12.0*z*(-x + xi) - 3.0*(x - xi)*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5))*(z**2 + (eta - y)**2 + (x - xi)**2)**22.5 + z*(eta - y)**2*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(6.0*z*(-x + xi) - 3.0*(x - xi)*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5))*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 - 3.0*z*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + 1.0*(-x + xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0 + 1.0*(-x + xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**24.5 + 3.0*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**23.0 + (x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**23.5*(2.0*z**2 + 1.0*z*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 3.0*(eta - y)**2) + (x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(1.0*z*(2.0*z + 2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5) + 2.0*(eta - y)**2)*(z**2 + (eta - y)**2 + (x - xi)**2)**24.0)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4
    d4Omega_dz2dxdy = lambda xi,eta: (z**2 + (eta - y)**2 + (x - xi)**2)**(-36.5)*(-15.0*z**2*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**33.0 + z**2*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(8.0*z + 6.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**33.5 + 6.0*z*(eta - y)*(x - xi)*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2 + (eta - y)**2 + (x - xi)**2)**34.5 + z*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(10.0*z + 2.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**34.0 + z*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(15.0*z**2*(z**2 + (eta - y)**2 + (x - xi)**2)**2.5 - 3.0*(z**2 + (eta - y)**2 + (x - xi)**2)**3.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**30.5 - 6.0*z*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**3*(z**2 + (eta - y)**2 + (x - xi)**2)**34.0 + 1.0*z*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(z**2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5 - (z**2 + (eta - y)**2 + (x - xi)**2)**1.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**33.0 + 2.0*z*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5 - (z**2 + (eta - y)**2 + (x - xi)**2)**1.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**33.5 + 1.0*z*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*((eta - y)*(3.0*z**2*(x - xi)*(z**2 + (eta - y)**2 + (x - xi)**2)**1.5 + 1.0*(-x + xi)*(z**2 + (eta - y)**2 + (x - xi)**2)**2.5) + (x - xi)*(3.0*z**2*(eta - y)*(z**2 + (eta - y)**2 + (x - xi)**2)**1.5 + 1.0*(-eta + y)*(z**2 + (eta - y)**2 + (x - xi)**2)**2.5))*(z**2 + (eta - y)**2 + (x - xi)**2)**32.0 - 4.0*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(1.0*z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**35.0 + 3.0*(eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4*(z**2 + (eta - y)**2 + (x - xi)**2)**34.0 - (eta - y)*(x - xi)*(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**2*(6.0*z + 2.0*(z**2 + (eta - y)**2 + (x - xi)**2)**0.5)*(z**2 + (eta - y)**2 + (x - xi)**2)**34.5)/(z + (z**2 + (eta - y)**2 + (x - xi)**2)**0.5)**4
    
    
    #Integrationen
    d3Qx_dx3 = simps(simps(q_x*d3Omega_dx3(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qx_dz3 = simps(simps(q_x*d3Omega_dz3(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qx_dx2dy = simps(simps(q_x*d3Omega_dx2dy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qx_dx2dz = simps(simps(q_x*d3Omega_dx2dz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qx_dy2dx = simps(simps(q_x*d3Omega_dy2dx(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qx_dz2dx = simps(simps(q_x*d3Omega_dz2dx(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qx_dz2dy = simps(simps(q_x*d3Omega_dz2dy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qx_dxdydz = simps(simps(q_x*d3Omega_dxdydz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    
    d3Qy_dy3 = simps(simps(q_y*d3Omega_dy3(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qy_dz3 = simps(simps(q_y*d3Omega_dz3(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qy_dx2dy = simps(simps(q_y*d3Omega_dx2dy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qy_dy2dx = simps(simps(q_y*d3Omega_dy2dx(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qy_dy2dz = simps(simps(q_y*d3Omega_dy2dz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qy_dz2dx = simps(simps(q_y*d3Omega_dz2dx(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qy_dz2dy = simps(simps(q_y*d3Omega_dz2dy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qy_dxdydz = simps(simps(q_y*d3Omega_dxdydz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    
    d3Qz_dz3 = simps(simps(q_z*d3Omega_dz3(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qz_dx2dz = simps(simps(q_z*d3Omega_dx2dz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qz_dy2dz = simps(simps(q_z*d3Omega_dy2dz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qz_dz2dx = simps(simps(q_z*d3Omega_dz2dx(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qz_dz2dy = simps(simps(q_z*d3Omega_dz2dy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d3Qz_dxdydz = simps(simps(q_z*d3Omega_dxdydz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    
    
    d4Qx_dx3dz = simps(simps(q_x*d4Omega_dx3dz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qx_dz3dx = simps(simps(q_x*d4Omega_dz3dx(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qx_dx2dz2 = simps(simps(q_x*d4Omega_dx2dz2(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qx_dx2dydz = simps(simps(q_x*d4Omega_dx2dydz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qx_dy2dxdz = simps(simps(q_x*d4Omega_dy2dxdz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qx_dz2dxdy = simps(simps(q_x*d4Omega_dz2dxdy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    
    d4Qy_dy3dz = simps(simps(q_y*d4Omega_dy3dz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qy_dz3dy = simps(simps(q_y*d4Omega_dz3dy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qy_dy2dz2 = simps(simps(q_y*d4Omega_dy2dz2(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qy_dx2dydz = simps(simps(q_y*d4Omega_dx2dydz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qy_dy2dxdz = simps(simps(q_y*d4Omega_dy2dxdz(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qy_dz2dxdy = simps(simps(q_y*d4Omega_dz2dxdy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    
    d4Qz_dz4 = simps(simps(q_z*d4Omega_dz4(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qz_dz3dx = simps(simps(q_z*d4Omega_dz3dx(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qz_dz3dy = simps(simps(q_z*d4Omega_dz3dy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qz_dx2dz2 = simps(simps(q_z*d4Omega_dx2dz2(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qz_dy2dz2 = simps(simps(q_z*d4Omega_dy2dz2(xi_vec,eta_vec), xi_linsp), eta_linsp)
    d4Qz_dz2dxdy = simps(simps(q_z*d4Omega_dz2dxdy(xi_vec,eta_vec), xi_linsp), eta_linsp)
    
    
    
    #Partielle Ableitungen der Komponenten von u; Quelle: Nikas, Boussinesq-Cerruti functions and a simple technique
    dux_dx = ((1+nu)/(2*np.pi*E))*(2*d3Qx_dz2dx - d3Qz_dx2dz + 2*nu*(d3Qx_dx3 + d3Qy_dx2dy + d3Qz_dx2dz) - z*(d4Qx_dx3dz + d4Qy_dx2dydz + d4Qz_dx2dz2))
    duy_dy = ((1+nu)/(2*np.pi*E))*(2*d3Qy_dz2dy - d3Qz_dy2dz + 2*nu*(d3Qx_dy2dx + d3Qy_dy3 + d3Qz_dy2dz) - z*(d4Qx_dy2dxdz + d4Qy_dy3dz + d4Qz_dy2dz2))
    duz_dz = ((1+nu)/(2*np.pi*E))*(d3Qz_dz3 - 2*nu*(d3Qx_dz2dx + d3Qy_dz2dy + d3Qz_dz3) - z*(d4Qx_dz3dx + d4Qy_dz3dy + d4Qz_dz4))
    
    dux_dy = ((1+nu)/(2*np.pi*E))*(2*d3Qx_dz2dy - d3Qz_dxdydz + 2*nu*(d3Qx_dx2dy + d3Qy_dy2dx + d3Qz_dxdydz) - z*(d4Qx_dx2dydz + d4Qy_dy2dxdz + d4Qz_dz2dxdy))
    dux_dz = ((1+nu)/(2*np.pi*E))*(2*d3Qx_dz3 - d3Qz_dz2dx + (2*nu-1)*(d3Qx_dx2dz + d3Qy_dxdydz + d3Qz_dz2dx) - z*(d4Qx_dx2dz2 + d4Qy_dz2dxdy + d4Qz_dz3dx))
    duy_dx = ((1+nu)/(2*np.pi*E))*(2*d3Qy_dz2dx - d3Qz_dxdydz + 2*nu*(d3Qx_dx2dy + d3Qy_dy2dx + d3Qz_dxdydz) - z*(d4Qx_dx2dydz + d4Qy_dy2dxdz + d4Qz_dz2dxdy))
    duy_dz = ((1+nu)/(2*np.pi*E))*(2*d3Qy_dz3 - d3Qz_dz2dy + (2*nu-1)*(d3Qx_dxdydz + d3Qy_dy2dz + d3Qz_dz2dy) - z*(d4Qx_dz2dxdy + d4Qy_dy2dz2 + d4Qz_dz3dy))
    duz_dx = ((1+nu)/(2*np.pi*E))*(d3Qz_dz2dx + (1-2*nu)*(d3Qx_dx2dz + d3Qy_dxdydz + d3Qz_dz2dx) - z*(d4Qx_dx2dz2 + d4Qy_dz2dxdy + d4Qz_dz3dx))
    duz_dy = ((1+nu)/(2*np.pi*E))*(d3Qz_dz2dy + (1-2*nu)*(d3Qx_dxdydz + d3Qy_dy2dz + d3Qz_dz2dy) - z*(d4Qx_dz2dxdy + d4Qy_dy2dz2 + d4Qz_dz3dy))
    
    
    
    #Normalspannungen
    sigma_x = (E/(1+nu))*((nu/(1-2*nu))*(dux_dx + duy_dy + duz_dz) + dux_dx)
    sigma_y = (E/(1+nu))*((nu/(1-2*nu))*(dux_dx + duy_dy + duz_dz) + duy_dy)
    sigma_z = (E/(1+nu))*((nu/(1-2*nu))*(dux_dx + duy_dy + duz_dz) + duz_dz)
    
    #Scherspannungen
    tau_xy = (E/(2+2*nu))*(dux_dy + duy_dx)
    tau_xz = (E/(2+2*nu))*(dux_dz + duz_dx)
    tau_yz = (E/(2+2*nu))*(duy_dz + duz_dy)
    
    #Vergleichsspannung GEH
    sigma_vGEH = ((1/2)*((sigma_x - sigma_y)**2 + (sigma_z - sigma_x)**2 + (sigma_y - sigma_z)**2 + 6*(tau_xy**2 + tau_xz**2 + tau_yz**2)))**(1/2)
    
    #Der in "vergleichsspannung_max" verwendete Optimierungsalgorithmus sucht das globale Minimum -> Das Vorzeichen muss vertauscht werden, um das Spannungsmaximum zu finden
    opt_vGEH = -sigma_vGEH
    
    # print()
    # print("vGEH =", -opt_vGEH)
    # print("x =", round(x,6), "\ty =", round(y,6), "\tz =", round(z,6))
    return opt_vGEH



def vergleichsspannung_max(E_1, nu_1, E_2, nu_2, a, b, q_z, q_x, q_y, xi_linsp, eta_linsp):
    """Dient zum finden des globalen Spannungsmaximums im Materialinneren"""
    
    #t_start, t_aktuell und t_delta dienen der Ausgabe der Rechendauer
    t_start = time.time()
    
    #Reduzierte Werkstoffkennwerte
    E = 1/((1/E_1)+(1/E_1))
    nu = (nu_1+nu_2)/2
    
    #Vektoren für die Integration
    xi_vec = xi_linsp[None, :]
    eta_vec = eta_linsp[:, None]
    
    #Randbedingungen für die Optimierung
    z_min = 1.5*(xi_linsp[1]-xi_linsp[0]) 
    bnds = ((-a, a), (-b, b), (z_min, a))  # Innerhalb dieses 3D-Raumes wird nach dem Spannungsmaximum gesucht
    
    "Suchen des globalen Spannungsmaximums mittels Optimierungsalgorithmus differential_evolution"
    #Aufruf vergleichsspannung
    solution = optimize.differential_evolution(vergleichsspannung, bnds, args=(E, nu, q_z, q_x, q_y, xi_linsp, eta_linsp, xi_vec, eta_vec), strategy='best1bin', maxiter=1000, popsize=9, tol=0.01, mutation=(0.4, 0.8), recombination=0.7, seed=None, callback=None, disp=False, polish=True, init='latinhypercube', atol=0, updating='immediate', workers=1) #, constraints=())
    
    #Erneute Vorzeichenänderung, um aus dem optimierten Minimum die maximale Vergleichsspannung zurückzugewinnen
    vGEH_max = -solution.fun
    
    #Position des Spannungsmaximums
    x_opt, y_opt, z_opt = solution.x
    
    #Ausgaben in der Konsole bezüglich Konvergenz und Ergebnis der Optimierung
    print()
    print(solution.message)
    t_aktuell = time.time()
    t_delta = t_aktuell - t_start
    print("Optimiert in", round(t_delta,1), "Sekunden")
    print("vGEH_max = %g" % (vGEH_max))
    print("x_opt = %g\ty_opt = %g\tz_opt = %g" % (x_opt, y_opt, z_opt))
    
    return vGEH_max, x_opt, y_opt, z_opt



def vergleichsspannung_plot(vGEH_max, x_opt, y_opt, z_opt, E_1, nu_1, E_2, nu_2, a, b, q_z, q_x, q_y, xi_linsp, eta_linsp):
    """Dient zum Berechnen der Spannungen im Materialinneren, welche später geplottet werden sollen"""
    
    #t_start, t_aktuell und t_delta dienen der Ausgabe der Rechendauer
    t_start = time.time()
    
    #Reduzierte Werkstoffkennwerte
    E = 1/((1/E_1)+(1/E_1))
    nu = (nu_1+nu_2)/2
    
    #Vektoren für die Integration
    xi_vec = xi_linsp[None, :]
    eta_vec = eta_linsp[:, None]
    
    #Angabe der für den Plot zu berechnenden Datenpunkte
    x_min_f = -1    # Faktor für die minimale x-Position
    x_max_f = 1     # Faktor für die maximale x-Position
    x_n = 20        # Anzahl der diskretisierten Abschnitte in x-Richtung
    
    y_min_f = -1    # y und z analog zu x
    y_max_f = 1
    y_n = 20
    
    z_min_f = 1.5*(xi_linsp[1]-xi_linsp[0])/a
    #z_min_f = 0.05
    z_max_f = 0.8
    z_n = 12
    
    #Erstellen der diskreten Datenpunkte in jeder Koordinatenrichtung anhand der zuvor festgelegten Werte
    x_linsp = np.linspace(x_min_f*a, x_max_f*a, x_n)
    y_linsp = np.linspace(y_min_f*b, y_max_f*b, y_n)
    z_linsp = np.linspace(z_min_f*a, z_max_f*a, z_n)
    
    
    "Berechnung der Vergleichsspannnungen für Plot in x-z-Ebene"
    #Initialisierung der Matrix, in die später die einzelnen vGEH geschrieben werden
    vGEH_zx = np.zeros((len(z_linsp),len(x_linsp)))
    
    #Berechnen der Vergleichsspannungen an den festgelegten Koordinaten durch den Aufruf der obigen Funktion "vergleichsspannung"
    for i in range(len(z_linsp)):
        for j in range(len(x_linsp)):
            koords = x_linsp[j], y_opt, z_linsp[i]
            vGEH_zx[i,j] = -vergleichsspannung(koords, E, nu, q_z, q_x, q_y, xi_linsp, eta_linsp, xi_vec, eta_vec)
    
    
    "Berechnung der Vergleichsspannnungen für Plot in y-z-Ebene"
    #Initialisierung der Matrix, in die später die einzelnen vGEH geschrieben werden
    vGEH_zy = np.zeros((len(z_linsp),len(y_linsp)))
    
    #Berechnen der Vergleichsspannungen an den festgelegten Koordinaten durch den Aufruf der obigen Funktion "vergleichsspannung"
    for i in range(len(z_linsp)):
        for k in range(len(y_linsp)):
            koords = x_opt, y_linsp[k], z_linsp[i]
            vGEH_zy[i,k] = -vergleichsspannung(koords, E, nu, q_z, q_x, q_y, xi_linsp, eta_linsp, xi_vec, eta_vec)
    
    # Ausgabe der Rechendauer in der Konsole
    t_aktuell = time.time()
    t_delta = t_aktuell - t_start
    print()
    print("Daten für Plots des Spannungsfeldes sind berechnet")
    print("Berechnet in", round(t_delta,1), "Sekunden")
    
    return x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy

