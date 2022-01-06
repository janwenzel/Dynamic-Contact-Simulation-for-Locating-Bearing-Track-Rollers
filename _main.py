# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:46:16 2020

@author: Christoph Bienefeld, Jan Wenzel, Alex Kretschmer

Dieses Skript beinhaltet die Hauptfunktion von Programm V1.

"""

import numpy as np
from scipy import optimize
import sys

#Selbst geschriebene Skripte importieren
import Hertz
import Kinematik
import Plots
import Potentialtheorie
import Spezialfall
import Woehlerlinie
import Materialskript


def optimierung_kraefte(r_delta0):
    """Übergeordnete Funktion, welche zum Finden des Kräftegleichgewichts genutzt wird"""
    
    #Festlegen der globalen Variablen
    global F_T1, F_T2, F_N1, F_N2, a1, b1, p_max1, a2, b2, p_max2, x_C_0spin1, x_C_0spin2, F_A1, M_R1, F_A2, M_R2, M_x_R_ist, M_x_R_delta, F_L, r_delta
    
    #for-Schleife liefert ausreichende Konvergenz für das Kräftegleichgewichtes aus Normal- und Tangentialkräften
    for ii in range(8):
        #Kräftegleichgewicht
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
        
        #Berechnen der x_C-Komponente des Bohrmittelpunktes für beide Kontakte
        x_C_0spin1 = r_delta0/np.sin(alpha_rad)
        x_C_0spin2 = -x_C_0spin1
        
        #Aufruf Kinematik.kontaktkraefte mit kurzer Rückgabe
        F_T1, F_A1, M_R1 = Kinematik.kontaktkraefte(n_x_C, a1, b1, p_max1, x_C_0spin1, alpha_rad1, mu, r_KonturFuehrung, r_Abroll, showPlots)
        F_T2, F_A2, M_R2 = Kinematik.kontaktkraefte(n_x_C, a2, b2, p_max2, x_C_0spin2, alpha_rad2, mu, r_KonturFuehrung, r_Abroll, showPlots)
        
        #Momentengleichgewicht in axialer Richtung
        M_x_R_ist = -(M_R1-M_R2)*np.sin(alpha_rad) - (F_A1+F_A2)*r_Abroll
        
        #Kräftegleichgewicht in Rollrichtung
        F_L = -F_A1-F_A2
        
        #r_delta ist globale Variable, sodass auf diese auch außerhalb der Funktion optimierung_kraefte zugegriffen werden kann
        r_delta = r_delta0
        
        #Mithilfe der Funktion optimierung_kraefte soll die Nullstelle von M_x_R_delta gefunden werden
        M_x_R_delta = M_x_R_ist - M_x_R_soll
        
    return M_x_R_delta




"Eingaben"

#Anzahl der diskreten Berechnungspunkte in x-Richtung
n_x_C = 40      #Werte zwischen 40 (minimal 20) und 80 (maximal 120) sind besonders empfehlenswert


#Rollenbelastung - Die zugrundeliegenden Kraftrichtungsdefinitionen können der Ausarbeitung (MT_Bienefeld) entnommen werden
#F_Motor > 0 -> gebremste Rolle
#F_Motor < 0 -> angetriebene Rolle

"--------NEU---------"
v_D = 60        #Durschnittliche Geschwindigkeit in [m/min]
"--------------------"     

F_Motor = 0            # Motorkraft [N]
F_ax = 0               # Axialkraft [N]        #positive Axialkraft zeigt in Rollrichtung rechts
F_rad = -210             # Radialkraft [N]      #negative Radialkraft sorgt für die Vorspannung der Rolle

mu = 0.3                # Reibungskoeffizient [-]


#Material der Laufrolle

"-----------NEU----------"        
Material_rolle="1.7225-unvergütet"

Material=Material_rolle

#Entnahme der Informationen über den Werkstoff auf "Materialskript"
E_Rolle,nu_Rolle,R_m_Rolle,sigma_D_Rolle,N_Rm_Rolle,N_k_Rolle=Materialskript.Auswahl(Material_rolle)


#Fuehrung
Material_Schiene="1.7225-unvergütet"  #!!! POM wurde im Materialskript editiert!!!

#Entnahme der Informationen über den Werkstoff auf "Materialskript"

E_Fuehrung,nu_Fuehrung,R_m_Schiene,sigma_D_Schiene,N_Rm_Schiene,N_k_Schiene=Materialskript.Auswahl(Material_Schiene)

"-------------"

#Geometrie des Gotik-Profils
alpha_deg = 30              # Kontaktwinkel [°]

A_T = 21.75                 # Abstand zwischen Rollenachse und Führungsachse [mm]
r_KonturRolle = -7          # Konturradius der Rolle in [mm]

r_KruemmungFuehrung = 39.5    # Kurven-Krümmungsradius der Führung [mm]
r_KonturFuehrung = 6        # Konturradius der Führung [mm]

Materialspannungsberechnung = "on"  # Wenn hier nicht "on" steht, wird die Berechnung der Vergleichsspannungen im Materialinneren ausgelassen


"Berechnung"

#Aufruf Hertz.parameter
E_ers, sumk, xi_IST, eta_IST, k_H, n_H, alpha_rad, r_Abroll, h_N, h_T = Hertz.parameter(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, alpha_deg, A_T, r_KonturRolle, r_KruemmungFuehrung, r_KonturFuehrung)

#alpha für die beiden Kontakte 1 (in Rollrichtung rechts) und 2 (in Rollrichtung links)
alpha_rad1 = alpha_rad
alpha_rad2 = -alpha_rad

#Maximale Last, bei der die Rolle vollständig durchdreht
F_L_max_abs = abs(mu*F_rad/np.cos(alpha_rad))


#Startbedingungen für optimierung_kraefte:
M_x_R_soll = F_Motor*r_Abroll   # Gefordertes Rollendrehmoment in axialer Richtung
F_T1 = 0                        # Die Tangentialkräfte an beiden Kontakten werden 
F_T2 = 0                        #   als Startbedingung zu null gesetzt
showPlots = False               # Parameter, welcher die Rückgabeart der Funktion Kinematik.kontaktkraefte bestimmt


try:        # Wenn die innerhalb von "try:" stehende Berechnung scheitert, dreht die Rolle durch -> es geht weiter bei "except ValueError:"
    print()
    print("Kräfte Optimierung:")
    #Aufruf optimierung_kraefte
    optimize.brentq(optimierung_kraefte, -r_Abroll*2/3, r_Abroll*2/3, maxiter=250)
    
    print("r_delta = ", r_delta)
    print("M_x_R_delta = ", M_x_R_delta)
    
    #Berechnung und Ausgabe von Schlupf und Reibmoment
    Schlupf = r_delta/r_Abroll
    Reibmoment = (M_R1-M_R2)*np.sin(alpha_rad)/1000     #Reibmoment bezüglich der Rollenachse in Nm
    print()
    print("Schlupf = ", round(Schlupf,8))
    print("Reibmoment bezüglich Rollenachse = ", round(Reibmoment,8), "Nm")
    print()
    
    #Ausgabe in der Konsole gibt Feedback zur Berechnung
    if M_x_R_soll != 0:
        M_x_R_fehler = abs(M_x_R_delta/M_x_R_soll)*100  #in Prozent
        if M_x_R_fehler > 5:
            print("--- Achtung: Berechnung der Kontaktkräfte ist ungenau! ---")
        else:
            print("Berechnung der Kontaktkräfte erfolgreich")
        print("Der Fehler beträgt", round(M_x_R_fehler,3),"%")
        print()
    
    
    
    "Plots"
    
    showPlots = True            # Parameter, welcher die Rückgabeart der Funktion Kinematik.kontaktkraefte bestimmt
    #Aufruf Kinematik.kontaktkraefte mit ausführlicher Rückgabe
    F_T1, F_A1, x_C_linsp1, y_C_linsp1, x_C_mesh1, y_C_mesh1, v_rel_x_C1, v_rel_y_C1, v_rel_res1, v_rel_x_C_dir1, v_rel_y_C_dir1, v_norm_x_C1, v_norm_y_C1, v_norm_res1, p1, tau_res1, tau_x_C1, tau_y_C1, dP2dxdy1, dPdy1, P1, dW2dxdy1, dWdy1, W1, M_R1, y_C_0spin1 = Kinematik.kontaktkraefte(n_x_C, a1, b1, p_max1, x_C_0spin1, alpha_rad1, mu, r_KonturFuehrung, r_Abroll, showPlots)
    F_T2, F_A2, x_C_linsp2, y_C_linsp2, x_C_mesh2, y_C_mesh2, v_rel_x_C2, v_rel_y_C2, v_rel_res2, v_rel_x_C_dir2, v_rel_y_C_dir2, v_norm_x_C2, v_norm_y_C2, v_norm_res2, p2, tau_res2, tau_x_C2, tau_y_C2, dP2dxdy2, dPdy2, P2, dW2dxdy2, dWdy2, W2, M_R2, y_C_0spin2 = Kinematik.kontaktkraefte(n_x_C, a2, b2, p_max2, x_C_0spin2, alpha_rad2, mu, r_KonturFuehrung, r_Abroll, showPlots)
    
    #Aufruf Plots.multiplot
    Plots.multiplot(F_N1, F_T1, F_A1, x_C_linsp1, y_C_linsp1, x_C_mesh1, y_C_mesh1, v_rel_x_C1, v_rel_y_C1, v_rel_res1, v_rel_x_C_dir1, v_rel_y_C_dir1, p1, tau_res1, tau_x_C1, tau_y_C1, dP2dxdy1, dPdy1, P1, dW2dxdy1, dWdy1, W1, F_N2, F_T2, F_A2, x_C_linsp2, y_C_linsp2, x_C_mesh2, y_C_mesh2, v_rel_x_C2, v_rel_y_C2, v_rel_res2, v_rel_x_C_dir2, v_rel_y_C_dir2, p2, tau_res2, tau_x_C2, tau_y_C2, dP2dxdy2, dPdy2, P2, dW2dxdy2, dWdy2, W2, a1, b1, a2, b2, x_C_0spin1, y_C_0spin1, x_C_0spin2, y_C_0spin2)
    
    #Aufruf Plots.print_kraefte
    Plots.print_kraefte(F_N1, F_T1, F_A1, F_N2, F_T2, F_A2, h_N, h_T)
    
    
    "-----------NEU-----------"
    # Berechnung der maximal ertragbaren Schwingspielzahl für die Laufrolle für Kontakt 1
    R_m=R_m_Rolle/3**(1/2)
    sigma_D=R_m_Rolle*0.3
    N_Rm=N_Rm_Rolle
    N_k=N_k_Rolle
    
    tau_max=p_max1*mu
    
    r_Abroll_ist=r_Abroll
    
    N,t_max=Woehlerlinie.SchwingspielzahlMAX_Schub(tau_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)

    N1_maxSchub_Rolle=round(N,2)
    t1_maxSchub_Rolle=round(t_max,2)

    # Berechnung der maximal ertragbaren Schwingspielzahl für die Schiene für Kontakt 1
    R_m=R_m_Schiene/3**(1/2)
    sigma_D=R_m_Schiene*0.3
    N_Rm=N_Rm_Schiene
    N_k=N_k_Schiene
    
    r_Abroll_ist=r_KruemmungFuehrung
    
    N,t_max=Woehlerlinie.SchwingspielzahlMAX_Schub(tau_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)              
    
    N1_maxSchub_Schiene=round(N,2)
    t1_maxSchub_Schiene=round(t_max,2)    
    
    
    # Berechnung der maximal ertragbaren Schwingspielzahl für die Laufrolle für Kontakt 2
    R_m=R_m_Rolle/3**(1/2)
    sigma_D=R_m_Rolle*0.3
    N_Rm=N_Rm_Rolle
    N_k=N_k_Rolle
    
    tau_max=p_max2*mu
    
    r_Abroll_ist=r_Abroll
    
    N,t_max=Woehlerlinie.SchwingspielzahlMAX_Schub(tau_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)

    N2_maxSchub_Rolle=round(N,2)
    t2_maxSchub_Rolle=round(t_max,2)

    # Berechnung der maximal ertragbaren Schwingspielzahl für die Schiene für Kontakt 2
    R_m=R_m_Schiene/3**(1/2)
    sigma_D=R_m_Schiene*0.3
    N_Rm=N_Rm_Schiene
    N_k=N_k_Schiene
    
    r_Abroll_ist=r_KruemmungFuehrung
    
    N,t_max=Woehlerlinie.SchwingspielzahlMAX_Schub(tau_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)              
    
    N2_maxSchub_Schiene=round(N,2)
    t2_maxSchub_Schiene=round(t_max,2)
    "-----------------------------------"
    
    
    "Potentialtheorie"
    
    if Materialspannungsberechnung == "on":
        print()
        print("Starte Berechnung der maximalen Vergleichsspannung mittels Potentialtheorie")
        
        print("Berechnung für Kontakt 1, da dieser maximal belastet ist")
        #Suche globales Spannungsmaximum: Aufruf Potentialtheorie.vergleichsspannung_max
        vGEH_max, x_opt, y_opt, z_opt = Potentialtheorie.vergleichsspannung_max(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a1, b1, p1, tau_x_C1, tau_y_C1, x_C_linsp1, y_C_linsp1)
        #Berechne die für die Plots benötigten Materialspannungen: Aufruf Potentialtheorie.vergleichsspannung_plot
        x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy = Potentialtheorie.vergleichsspannung_plot(vGEH_max, x_opt, y_opt, z_opt, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a1, b1, p1, tau_x_C1, tau_y_C1, x_C_linsp1, y_C_linsp1)
        
        "-----------NEU-----------"
        # Berechnung der maximal ertragbaren Schwingspielzahl für die Laufrolle für Kontakt 1
        R_m=R_m_Rolle
        sigma_D=sigma_D_Rolle
        N_Rm=N_Rm_Rolle
        N_k=N_k_Rolle
        
        r_Abroll_ist=r_Abroll
        
        N,t_max=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)

        N1_max_Rolle=round(N,2)
        t1_max_Rolle=round(t_max,2)

        # Berechnung der maximal ertragbaren Schwingspielzahl für die Schiene für Kontakt 1
        R_m=R_m_Schiene
        sigma_D=sigma_D_Schiene
        N_Rm=N_Rm_Schiene
        N_k=N_k_Schiene
        
        r_Abroll_ist=r_KruemmungFuehrung
        
        N,t_max=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)              
        
        N1_max_Schiene=round(N,2)
        t1_max_Schiene=round(t_max,2)

        "------------------------"
        #Aufruf Plots.plot_vergleichsspannung
        Plots.plot_vergleichsspannung(x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy, x_opt, y_opt, z_opt)
        
        print("Berechnung für Kontakt 2, da dieser maximal belastet ist")
        #Suche globales Spannungsmaximum: Aufruf Potentialtheorie.vergleichsspannung_max
        vGEH_max, x_opt, y_opt, z_opt = Potentialtheorie.vergleichsspannung_max(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a2, b2, p2, tau_x_C2, tau_y_C2, x_C_linsp2, y_C_linsp2)
        #Berechne die für die Plots benötigten Materialspannungen: Aufruf Potentialtheorie.vergleichsspannung_plot
        x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy = Potentialtheorie.vergleichsspannung_plot(vGEH_max, x_opt, y_opt, z_opt, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a2, b2, p2, tau_x_C2, tau_y_C2, x_C_linsp2, y_C_linsp2)
         
        "-----------NEU-----------"
        # Berechnung der maximal ertragbaren Schwingspielzahl für die Laufrolle für Kontakt 1
        R_m=R_m_Rolle
        sigma_D=sigma_D_Rolle
        N_Rm=N_Rm_Rolle
        N_k=N_k_Rolle
        
        r_Abroll_ist=r_Abroll
        
        N,t_max=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)

        N2_max_Rolle=round(N,2)
        t2_max_Rolle=round(t_max,2)

        # Berechnung der maximal ertragbaren Schwingspielzahl für die Schiene für Kontakt 1
        R_m=R_m_Schiene
        sigma_D=sigma_D_Schiene
        N_Rm=N_Rm_Schiene
        N_k=N_k_Schiene
        
        r_Abroll_ist=r_KruemmungFuehrung
        
        N,t_max=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist)              
        
        N2_max_Schiene=round(N,2)
        t2_max_Schiene=round(t_max,2)
        
        "------------------------"
        
        #Aufruf Plots.plot_vergleichsspannung
        Plots.plot_vergleichsspannung(x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy, x_opt, y_opt, z_opt)
    

#Wenn die Rolle durchdreht, scheitert obige Berechnung bei "try:" mit einem "ValueError" -> Es geht hier weiter
except ValueError:
    print("--- Achtung: Rolle dreht durch! ---")
    print("--- Programm wird mit Spezialfall.durchdrehen fortgesetzt ---")
    #Aufruf Spezialfall.durchdrehen
    Spezialfall.durchdrehen(n_x_C, F_Motor, F_ax, F_rad, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, mu, alpha_deg, A_T, r_KonturRolle, r_KruemmungFuehrung, r_KonturFuehrung, Materialspannungsberechnung)


"-------Ausgabe------"
#%%
N_Rolle_max_arr=np.array([N1_max_Rolle,N2_max_Rolle])
N_Rolle_max=np.amin(N_Rolle_max_arr)

t_Rolle_max_arr=np.array([t1_max_Rolle,t2_max_Rolle])
t_Rolle_max=np.amin(t_Rolle_max_arr)

# Gleiches für Schub:
N_Rolle_maxSchub_arr=np.array([N1_maxSchub_Rolle,N2_maxSchub_Rolle])
N_Rolle_maxSchub=np.amin(N_Rolle_maxSchub_arr)

t_Rolle_maxSchub_arr=np.array([t1_maxSchub_Rolle,t2_maxSchub_Rolle])
t_Rolle_maxSchub=np.amin(t_Rolle_maxSchub_arr)

print("----------------------")
print("Durchschnittliche Geschwindigkeit = ",v_D," m/min")
print("----------------------")

print("Überrollungen Laufrolle = ",N_Rolle_max)
print("Dauer bis Anriss Laufrolle = ",t_Rolle_max, " Stunden")

print("Überrollungen Laufrolle Schub = ",N_Rolle_maxSchub)
print("Dauer bis Oberflächenschäden Laufrolle Schub = ",t_Rolle_maxSchub, " Stunden")
print("----------------------")

N_Schiene_max_arr=np.array([N1_max_Schiene,N2_max_Schiene])
N_Schiene_max=np.amin(N_Schiene_max_arr)

t_Schiene_max_arr=np.array([t1_max_Schiene,t2_max_Schiene])
t_Schiene_max=np.amin(t_Schiene_max_arr)

# Gleiches für Schub:
N_Schiene_maxSchub_arr=np.array([N1_maxSchub_Schiene,N2_maxSchub_Schiene])
N_Schiene_maxSchub=np.amin(N_Schiene_maxSchub_arr)

t_Schiene_maxSchub_arr=np.array([t1_maxSchub_Schiene,t2_maxSchub_Schiene])
t_Schiene_maxSchub=np.amin(t_Schiene_maxSchub_arr)

print("Überrollungen Schiene = ",N_Schiene_max)
print("Dauer bis Anriss Schiene = ",t_Schiene_max, " Stunden")

print("Überrollungen Schiene Schub = ",N_Schiene_maxSchub)
print("Dauer bis Oberflächenschäden Schiene Schub = ",t_Schiene_maxSchub, " Stunden")     
        

