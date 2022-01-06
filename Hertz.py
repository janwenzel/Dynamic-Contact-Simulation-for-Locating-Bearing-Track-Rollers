# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:28:50 2020

@author: Jan Wenzel, edited by Christoph Bienefeld
"""

import numpy as np
from scipy.interpolate import interp1d


def pressung(F_N, sumk, E_ers, xi_IST, eta_IST):
    """Berechnet die beiden Halbachsen der Kontaktellipse a und b sowie die maximale Normalspannung an der Oberfläche p_max"""
    
    #Berechnung der Halbachsen
    a=xi_IST*((3*E_ers)/sumk)**(1/3)*F_N**(1/3)
    
    b=eta_IST*((3*E_ers)/sumk)**(1/3)*F_N**(1/3)
    
    #Berechnung der maximalen Pressung
    p_max=3/2*F_N/(np.pi*a*b)
    
    return a, b, p_max



def parameter(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, alpha_deg, A_T, r_KonturRolle, r_KruemmungFuehrung, r_KonturFuehrung):
    """Vorberechnungen von Geometrie-Größen und Hertz-Parametern"""
    
    #Ersatz-E-Modul
    E_ers=0.5*(((1-nu_Rolle**2)/E_Rolle)+((1-nu_Fuehrung**2)/E_Fuehrung)) #E_ers ersetzt (1-nu^2)/E
    
    #Geometrie-Daten:
    alpha_rad = (alpha_deg/180)*np.pi
    r_Abroll = A_T - r_KonturFuehrung*np.cos(alpha_rad)
    r_AbrollIst = r_Abroll/np.cos(alpha_rad)
    r_KurvenKruemmung = r_KruemmungFuehrung/np.cos(alpha_rad)
    
    #Hebelarme für die Berechnung des Kippmoments
    h_N = A_T*np.sin(alpha_rad)
    h_T = A_T*np.cos(alpha_rad) - r_KonturFuehrung
    
    
    #Die Radien R_1 und R_2 werden in Anlehnung an Hamrock, Dawson eingeführt.
    #Diese sind definiert als: 1/R_x=1/r_ax+1/r_bx und 1/R_y=1/r_ay+1/r_by
    #Hierbei steht a und b für die beiden Körper und x und y für die beiden
    #betrachteten Ebenen, in denen die Krümmungsradien liegen:
    R_1=1/(1/r_AbrollIst+1/r_KurvenKruemmung)
    R_2=1/(1/r_KonturRolle+1/r_KonturFuehrung)
    
    
    #----Festlegen der Indizierung der Krümmungsradien---
    #Der spannende Teil in diesem Skript! np.where ist quasi eine if, else Bedingung, die aber mit ...
    #...linspace kombiniert werden kann   
    r_11=np.where(1/R_1>1/R_2,r_AbrollIst, r_KonturRolle)
    r_21=np.where(1/R_1>1/R_2,r_KurvenKruemmung, r_KonturFuehrung)
    r_12=np.where(1/R_1>1/R_2,r_KonturRolle, r_AbrollIst)
    r_22=np.where(1/R_1>1/R_2, r_KonturFuehrung, r_KurvenKruemmung)
    
    
    #Berechnung der Krümmungswerte:
    k_11=1/r_11
    k_12=1/r_12
    k_21=1/r_21
    k_22=1/r_22
    
    
    sumk=k_11+k_12+k_21+k_22
    
    costau=(abs(k_11-k_12+k_21-k_22)/sumk)
    
    
    #Beiwerte Hertz
    
    #Eingabe der Tabellen - Werte entnommen aus Schaeffler Technisches Taschenbuch, S.255-256
    costau_table=  [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.9995]
    xi_table=      [1.07, 1.15, 1.24, 1.35, 1.48, 1.66, 1.91, 2.3, 3.09, 4.12, 7.76, 10.15, 23.95]
    eta_table=     [0.938, 0.879, 0.824, 0.771, 0.718, 0.664, 0.607, 0.544, 0.461, 0.396, 0.287, 0.251, 0.163]
    psiDIVxi_table=[0.997, 1.01, 1.02, 0.962, 0.938, 0.904, 0.859, 0.792, 0.680, 0.577, 0.384, 0.320, 0.171]
    
    #Umwandeln der eingegebenen Liste als Array mit NumPy
    costau_array=np.array(costau_table)
    xi_array=np.array(xi_table)
    eta_array=np.array(eta_table)
    psiDIVxi_array=np.array(psiDIVxi_table)
    
    #Interpolieren der Tabellenwerte und erstellen der Funktionen
    funxi=interp1d(costau_array, xi_array)
    funeta=interp1d(costau_array, eta_array)
    funpsiDIVxi=interp1d(costau_array, psiDIVxi_array)
    
    #Berechnung der IST-Werte mit den interp. Funktionen
    xi_IST=funxi(costau)
    eta_IST=funeta(costau)
    psiDIVxi_IST=funpsiDIVxi(costau)
    
    #Berechnung der Steifigkeit
    k_H=1/(psiDIVxi_IST**(3/2)*((9*sumk*E_ers**2)/8)**(1/2))
    
    #Berechnung der Nachgiebigkeit
    n_H=psiDIVxi_IST*((9*sumk*E_ers**2)/8)**(1/3)
    
    
    return E_ers, sumk, xi_IST, eta_IST, k_H, n_H, alpha_rad, r_Abroll, h_N, h_T

