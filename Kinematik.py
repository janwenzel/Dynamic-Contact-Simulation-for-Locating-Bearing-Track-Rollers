# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 17:19:46 2020

@author: Christoph Bienefeld
"""

import numpy as np
from scipy.integrate import simps #Alternativ: trapz


def kontaktkraefte(n_x_C, a, b, p_max, x_C_0spin, alpha_rad, mu, r_KonturFuehrung, r_Abroll, showPlots):
    """Berechnung der Spannungsverteilung an der Oberfläche sowie der Kontaktkräfte bei regulärem Rollen"""
    
    #Angabe der Rollenumdrehungen pro Minute, der Wert 60 entspricht einer Umdrehung pro Sekunde
    #Dadurch ensprechen alle Leistungen im Betrag den Energien pro Umdrehung
    drehzahl_Rolle = 60     # [U/min]
    
    
    "Berechnung"
    #Winkelgeschwindigkeit in rad/s umrechnen
    omega_Rolle = drehzahl_Rolle*2*np.pi/60
    
    #Rechteck um die Kontaktellipse diskretisieren
    x_C_linsp = np.linspace(-a, a, n_x_C)
    y_C_linsp = np.linspace(-b, b, int(round(n_x_C*b/a)))
    x_C_mesh, y_C_mesh = np.meshgrid(x_C_linsp, y_C_linsp)
    
    
    #y_C-Komponente des Bohrmittelpunktes berechnen
    x_C_spinParam = x_C_0spin/a
    if x_C_spinParam > -1 and x_C_spinParam < 1:
        y_C_spinParam = -(1 - (x_C_spinParam)**2)**(1/2) + 0
    else:
        y_C_spinParam = 0
    y_C_0spin = y_C_spinParam*b
    
    
    #y-Komponente des Bohrmittelpunktes auf das diskretisierte Rechteck abbilden
    y_C_delta_linsp = np.zeros(n_x_C)
    for ee in range(n_x_C):
        if x_C_linsp[ee] < x_C_0spin:
            y_C_delta_linsp[ee] = (x_C_linsp[ee] + a)*(y_C_0spin)/(x_C_0spin + a)
        else:
            y_C_delta_linsp[ee] = (x_C_linsp[ee] - a)*(y_C_0spin)/(x_C_0spin - a)
    y_C_delta_linsp = np.nan_to_num(y_C_delta_linsp)
    y_C_delta_mesh = np.meshgrid(y_C_delta_linsp, y_C_linsp)[0]
    # y_C_delta_mesh = y_C_0spin    #Wenn diese Zeile aktiv ist, wird Annahme 2 deaktiviert (siehe Ausarbeitung der Masterarbeit)
    
    
    #Relativgeschwindigkeiten im Kontakt berechnen
    v_rel_x_C = np.sin(alpha_rad)*omega_Rolle*(y_C_mesh - y_C_delta_mesh)
    v_rel_y_C = omega_Rolle*((r_KonturFuehrung**2 - (x_C_mesh*np.cos(alpha_rad) + r_KonturFuehrung*np.sin(alpha_rad))**2)**(1/2) - r_KonturFuehrung*np.cos(alpha_rad) + x_C_0spin*np.sin(alpha_rad))
    v_rel_res = np.hypot(v_rel_x_C, v_rel_y_C)
        
    #Richtungen der Relativgeschwindigkeiten durch Normierung berechnen
    v_rel_x_C_dir = v_rel_x_C/abs(v_rel_res)
    v_rel_y_C_dir = v_rel_y_C/abs(v_rel_res)
    
    #Auf die Rollgeschwindigkeit normierte Relativgeschwindigkeiten
    v_norm_x_C = v_rel_x_C/(omega_Rolle*r_Abroll)
    v_norm_y_C = v_rel_y_C/(omega_Rolle*r_Abroll)
    v_norm_res = v_rel_res/(omega_Rolle*r_Abroll)
    
    #Pressungsverteilung im Kontakt ermitteln
    p = p_max*(1-(x_C_mesh/a)**2-(y_C_mesh/b)**2)**(1/2)
    p = np.nan_to_num(p)
    
    #Tangentialspannungen aus Reibwert, Pressung und Geschwindigkeitsrichtungen bestimmen
    tau_res = mu*p
    tau_x_C = -tau_res*v_rel_x_C_dir
    tau_y_C = -tau_res*v_rel_y_C_dir
    
    #Tangentialkraft und Antriebskraft durch Integration der Tangentialspannungen bestimmen
    T_x_C = simps(simps(tau_x_C, x_C_linsp), y_C_linsp)
    T_y_C = simps(simps(tau_y_C, x_C_linsp), y_C_linsp)
    
    
    #Reibleistung aus Tangentialspannungen und Relativgeschwindigkeiten ermitteln
    dP2dxdy = tau_res*v_rel_res
    dPdy = simps(dP2dxdy, y_C_linsp, axis = 0)
    P = simps(simps(dP2dxdy, x_C_linsp), y_C_linsp)
    
    #Reibenergien pro Rollweg
    dW2dxdy = tau_res*v_norm_res
    dWdy = simps(dW2dxdy, y_C_linsp, axis = 0)
    W = simps(simps(dW2dxdy, x_C_linsp), y_C_linsp)
    
    #Reibmoment berechnen
    d2M_z_Cdxdy = -tau_x_C*y_C_mesh + tau_y_C*x_C_mesh
    M_z_C = simps(simps(d2M_z_Cdxdy, x_C_linsp), y_C_linsp)
    
    "Festlegen der Rückgabevariablen: Entweder nur für die Optimierung oder für das Erzeugen der Plots"
    if showPlots:
        return T_x_C, T_y_C, x_C_linsp, y_C_linsp, x_C_mesh, y_C_mesh, v_rel_x_C, v_rel_y_C, v_rel_res, v_rel_x_C_dir, v_rel_y_C_dir, v_norm_x_C, v_norm_y_C, v_norm_res, p, tau_res, tau_x_C, tau_y_C, dP2dxdy, dPdy, P, dW2dxdy, dWdy, W, M_z_C, y_C_0spin
        
    else:
        return T_x_C, T_y_C, M_z_C
    


def kin_durchdrehen(F_L_sign, n_x_C, a, b, p_max, alpha_rad, mu, r_KonturFuehrung, r_Abroll):
    """Berechnung der Spannungsverteilung an der Oberfläche sowie der Kontaktkräfte bei durchdrehender Rolle"""
    
    #Rechteck um die Kontaktellipse diskretisieren
    x_C_linsp = np.linspace(-a, a, n_x_C)
    y_C_linsp = np.linspace(-b, b, int(round(n_x_C*b/a)))
    x_C_mesh, y_C_mesh = np.meshgrid(x_C_linsp, y_C_linsp)
    
    #Pressungsverteilung im Kontakt ermitteln
    p = p_max*(1-(x_C_mesh/a)**2-(y_C_mesh/b)**2)**(1/2)
    p = np.nan_to_num(p)
    
    #Tangentialspannungen aus Reibwert, Pressung und Geschwindigkeitsrichtungen bestimmen
    tau_res = -mu*p*F_L_sign
    tau_x_C = 0
    tau_y_C = -mu*p*F_L_sign
    
    #Tangentialkraft in Antriebsrichtung bestimmen
    T_y_C = simps(simps(tau_y_C, x_C_linsp), y_C_linsp)
    
    return T_y_C, x_C_linsp, y_C_linsp, x_C_mesh, y_C_mesh, p, tau_res, tau_x_C, tau_y_C

