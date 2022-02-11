
"""

Technical Universitiy Darmstadt, Germany
Institute for Product Development and Machine Elements

Defines input parameters and starts the simulation

The simulation will start when running this file.

To have access to the full ouput, you can start the "_main.py" script directly.

"""

def parameter():

    "Forces"
    #Rollenbelastung - Die zugrundeliegenden Kraftrichtungsdefinitionen können der Dokumentation entnommen werden
    #F_Motor > 0 -> Braking 
    #F_Motor < 0 -> Driving
    F_Motor = 30            # Driving Force [N]
    F_ax = 30               # Axial Force [N]        #positive Axialkraft zeigt in Rollrichtung rechts
    F_rad = -200             # Radial Force [N]      #negative Radialkraft sorgt für die Vorspannung der Rolle
    
    "Material"
    #Material der Laufrolle   
    Material_rolle="1.7225-unvergütet"
    #Material der Fuehrung
    Material_Schiene="1.7225-unvergütet"
    # Coefficient of Friction [-]
    mu = 0.3                
    
    "Geometry"    
    #Geometrie des Gotik-Profils
    alpha_deg = 30              # Kontaktwinkel [°]
    
    A_T = 21.75                 # Abstand zwischen Rollenachse und Führungsachse [mm]
    r_KonturRolle = -7          # Konturradius der Rolle in [mm]
    
    r_KruemmungFuehrung = 39.5    # Kurven-Krümmungsradius der Führung [mm]
    r_KonturFuehrung = 6        # Konturradius der Führung [mm]
    
    "Velocity"
    #Mean velocity [m/min] (For life-time prediction)
    v_D = 60 
    
    "Simulation Settings"
    #Anzahl der diskreten Berechnungspunkte in x-Richtung
    n_x_C = 40     #Werte zwischen 40 (minimal 20) und 80 (maximal 120) sind besonders empfehlenswert 
    #Auswahl-Optionen
    Materialspannungsberechnung = "on"  # Wenn hier nicht "on" steht, wird die Berechnung der Vergleichsspannungen im Materialinneren ausgelassen

    return n_x_C,v_D,F_Motor,F_ax,F_rad,mu,Material_rolle,Material_Schiene,alpha_deg,A_T,r_KonturRolle,r_KruemmungFuehrung,r_KonturFuehrung,Materialspannungsberechnung
            
#Runs the simulation
def run():
    import _main
    _main

run()