"""

Technical Universitiy Darmstadt, Germany
Institute for Product Development and Machine Elements

Dynamic Contact Simulation for Locating Bearing Track Rollers

Woehler Curve

"""

import numpy as np

def SchwingspielzahlMAX (vGEH_max,R_m,sigma_D,N_Rm,N_k,v_D,r_Abroll_ist):

    #Berechnung der Steigung
    k=np.log(N_Rm/N_k)/np.log(sigma_D/R_m)    
    
    #Berechnung der ertragbaren Schwingspielzahls 
    N=(N_k*(vGEH_max/sigma_D)**(-k))

    v_D=v_D/60*1000     #Umrechnung der Geschwindigkeit von m/min in mm/s    

    n=v_D*30/(np.pi*r_Abroll_ist) #Umdrehung in 1/min

    t_max=N/n/60      #Umrechnung der Zeit von min in h
         
            
    return (N,t_max)


