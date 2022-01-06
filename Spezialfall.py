# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:17:47 2020

@author: Christoph Bienefeld
"""

import numpy as np
import sys

import Hertz
import Kinematik
import Plots
import Potentialtheorie


def durchdrehen(n_x_C, F_Motor, F_ax, F_rad, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, mu, alpha_deg, A_T, r_KonturRolle, r_KruemmungFuehrung, r_KonturFuehrung, Materialspannungsberechnung):
    """Übergeordnete Funktion zur Berechnung der Spannungszustände bei durchdrehener Rolle"""
    
    #Aufruf Hertz.parameter
    E_ers, sumk, xi_IST, eta_IST, k_H, n_H, alpha_rad, r_Abroll, h_N, h_T = Hertz.parameter(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, alpha_deg, A_T, r_KonturRolle, r_KruemmungFuehrung, r_KonturFuehrung)
    
    #alpha für die beiden Kontakte 1 (in Rollrichtung rechts) und 2 (in Rollrichtung links)
    alpha_rad1 = alpha_rad
    alpha_rad2 = -alpha_rad
    
    #Vorzeichen der Rollenbelastung
    F_L_sign = np.sign(F_Motor)
    
    
    #Tangentialkräfte bei durchdrehender Rolle
    F_T1 = 0
    F_T2 = 0
    
    #Berechnung der Normalkräfte
    F_N1 = 0.5*((-F_ax - np.cos(alpha_rad)*(F_T1+F_T2))/(np.sin(alpha_rad)) + (-F_rad + np.sin(alpha_rad)*(F_T1-F_T2))/(np.cos(alpha_rad)))
    F_N2 = 0.5*((F_ax + np.cos(alpha_rad)*(F_T1+F_T2))/(np.sin(alpha_rad)) + (-F_rad + np.sin(alpha_rad)*(F_T1-F_T2))/(np.cos(alpha_rad)))
    
    #Abbruch bei Abheben des Rades
    if F_N1 < 0 or F_N2 < 0:
        print("--- Achtung: Rolle hebt ab! ---")
        print("--- Programm wird gestoppt ---")
        print()
        sys.exit(1)
    
    
    #Aufruf Hertz.pressung
    a1, b1, p_max1 = Hertz.pressung(F_N1, sumk, E_ers, xi_IST, eta_IST)
    a2, b2, p_max2 = Hertz.pressung(F_N2, sumk, E_ers, xi_IST, eta_IST)
    
    
    #Aufruf Kinematik.kin_durchdrehen
    F_A1, x_C_linsp1, y_C_linsp1, x_C_mesh1, y_C_mesh1, p1, tau_res1, tau_x_C1, tau_y_C1 = Kinematik.kin_durchdrehen(F_L_sign, n_x_C, a1, b1, p_max1, alpha_rad1, mu, r_KonturFuehrung, r_Abroll)
    F_A2, x_C_linsp2, y_C_linsp2, x_C_mesh2, y_C_mesh2, p2, tau_res2, tau_x_C2, tau_y_C2 = Kinematik.kin_durchdrehen(F_L_sign, n_x_C, a2, b2, p_max2, alpha_rad2, mu, r_KonturFuehrung, r_Abroll)
    
    #Momentengleichgewicht bei durchdrehender Rolle
    M_x_R = -(F_A1+F_A2)*r_Abroll
    F_Motor = M_x_R/r_Abroll
    print("Motorkraft F_Motor =", F_Motor, "N" )
    
    #Aufruf Plots.print_kraefte
    Plots.print_kraefte(F_N1, F_T1, F_A1, F_N2, F_T2, F_A2, h_N, h_T)
    
    
    "Potentialtheorie"
    
    if Materialspannungsberechnung == "on":
        #Aufruf Potentialtheorie
        print()
        print("Starte Berechnung der maximalen Vergleichsspannung mittels Potentialtheorie")
        print("Berechnung für Kontakt 1")#, da dieser maximal belastet ist")
        vGEH_max1, x_opt1, y_opt1, z_opt1 = Potentialtheorie.vergleichsspannung_max(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a1, b1, p1, tau_x_C1, tau_y_C1, x_C_linsp1, y_C_linsp1)
        #x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy = Potentialtheorie.vergleichsspannung_plot(vGEH_max, x_opt, y_opt, z_opt, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a1, b1, p1, tau_x_C1, tau_y_C1, x_C_linsp1, y_C_linsp1)
        print("Berechnung für Kontakt 2")#", da dieser maximal belastet ist")
        vGEH_max2, x_opt2, y_opt2, z_opt2 = Potentialtheorie.vergleichsspannung_max(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a2, b2, p2, tau_x_C2, tau_y_C2, x_C_linsp2, y_C_linsp2)
        #x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy = Potentialtheorie.vergleichsspannung_plot(vGEH_max, x_opt, y_opt, z_opt, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a2, b2, p2, tau_x_C2, tau_y_C2, x_C_linsp2, y_C_linsp2)
        #Plots.plot_vergleichsspannung(x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy, x_opt, y_opt, z_opt)
    
    # noch return hinzufügen?
    if Materialspannungsberechnung == "on":
        return F_N1, F_N2, F_T1, F_T2, F_A1, F_A2, p_max1, p_max2, tau_x_C1, tau_x_C2, tau_y_C1, tau_y_C2, a1, a2, b1, b2, vGEH_max1, x_opt1, y_opt1, z_opt1, vGEH_max2, x_opt2, y_opt2, z_opt2  
    else: 
        return F_N1, F_N2, F_T1, F_T2, F_A1, F_A2, p_max1, p_max2, tau_x_C1, tau_x_C2, tau_y_C1, tau_y_C2, a1, a2, b1, b2