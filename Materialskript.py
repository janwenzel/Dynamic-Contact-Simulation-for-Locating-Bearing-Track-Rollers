"""

Technical Universitiy Darmstadt, Germany
Institute for Product Development and Machine Elements

Dynamic Contact Simulation for Locating Bearing Track Rollers

Material Database

"""

# Auswahl: Hier werden alle Materialien mit ihren Eigenschaften aufgelistet
# Reibkoeffizient: Hier werden die jeweiligen Reibkoeffizienten für die Materialpaarung definiert


def Auswahl (Material):
    
    "----Stähle------"    
    
    if Material=="1.7225-unvergütet":
        
        E=210000
        nu=0.33
        
        R_m=720
        sigma_D=0.5*R_m
        N_Rm=10**4
        N_k=5*10**6        

    if Material=="1.7225-vergütet":
        
        E=210000
        nu=0.33
        
        R_m=1200
        sigma_D=0.5*R_m
        N_Rm=10**4
        N_k=5*10**6
    
    
    "----Aluminiumlegierungen-----"
    
    if Material=="ENAW7075":
        
        E=70000
        nu=0.33
        
        R_m=540
        sigma_D=0.5*R_m
        N_Rm=10**4
        N_k=5*10**6
    
    
    
    "----Kunststoffe-----"
    if Material=="PA6":
        
        E=5000
        nu=0.3
        
        R_m=60
        sigma_D=0.3*R_m
        N_Rm=10**4
        N_k=5*10**6


    if Material=="POM":
        
        E=3000
        nu=0.3
        
        R_m=80
        sigma_D=0.4*R_m
        N_Rm=10**4
        N_k=5*10**6
        
        
    
    return (E,nu,R_m,sigma_D,N_Rm,N_k)


# Werte kommen aus den genannten Quellen in "Quellen Reibkoeffizienten.txt"
def Reibkoeffizient (Material_Rolle, Material_Schiene, Reibkoeff_str):
    
    if Reibkoeff_str == "Realwert":
        mu_Metall_Metall = 0.19 # Popov
        mu_Metall_PA6 = 0.22 # triboplast.de
        mu_Metall_POM = 0.14 # triboplast.de
        mu_PA6_PA6 = 0.2 # Platzhalter
        mu_PA6_POM = 0.2 # Platzhalter
        mu_POM_POM = 0.2 # Platzhalter
    if Reibkoeff_str == "Oberwert":
        mu_Metall_Metall = 0.3 # schweizer-fn.de
        mu_Metall_PA6 = 0.3 # schweizer-fn.de
        mu_Metall_POM = 0.3 # Platzhalter
        mu_PA6_PA6 = 0.3 # Platzhalter
        mu_PA6_POM = 0.3 # Platzhalter
        mu_POM_POM = 0.3 # Platzhalter
    if Reibkoeff_str == "Unterwert":
        mu_Metall_Metall = 0.16 # Popov
        mu_Metall_PA6 = 0.1 # triboplast.de
        mu_Metall_POM = 0.06 # triboplast.de
        mu_PA6_PA6 = 0.1 # Platzhalter
        mu_PA6_POM = 0.1 # Platzhalter
        mu_POM_POM = 0.06 # Platzhalter
    
    if Material_Rolle == "1.7225-unvergütet" or Material_Rolle == "1.7225-vergütet" or Material_Rolle == "ENAW7075":
        if Material_Schiene == "1.7225-unvergütet" or Material_Schiene == "1.7225-vergütet" or Material_Schiene == "ENAW7075":
            mu = mu_Metall_Metall
        if Material_Schiene == "PA6":
            mu = mu_Metall_PA6
        if Material_Schiene == "POM":
            mu = mu_Metall_POM
    
    if Material_Rolle == "PA6":
        if Material_Schiene == "1.7225-unvergütet" or Material_Schiene == "1.7225-vergütet" or Material_Schiene == "ENAW7075":
            mu = mu_Metall_PA6
        if Material_Schiene == "PA6": 
            mu = mu_PA6_PA6
        if Material_Schiene == "POM":
            mu = mu_PA6_POM
    
    if Material_Rolle == "POM":
        if Material_Schiene == "1.7225-unvergütet" or Material_Schiene == "1.7225-vergütet" or Material_Schiene == "ENAW7075":
            mu = mu_Metall_POM
        if Material_Schiene == "PA6": 
            mu = mu_PA6_POM
        if Material_Schiene == "POM":
            mu = mu_POM_POM    
    
    return mu
    
    
    