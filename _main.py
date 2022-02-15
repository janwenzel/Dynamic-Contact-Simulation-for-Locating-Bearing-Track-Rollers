
"""

Technical Universitiy Darmstadt, Germany
Institute for Product Development and Machine Elements

Main-Script

"""

import numpy as np
from scipy import optimize
import sys

import Hertz
import Kinematik
import Plots
import Potentialtheorie
import Spezialfall
import Woehlerlinie
import Materialskript
from Input_and_Run_Sim import parameter


def optimierung_kraefte(r_delta0):
    """Overall function which is used to find the force equilibrium"""
    
    #Setting the global variables
    global F_T1, F_T2, F_N1, F_N2, a1, b1, p_max1, a2, b2, p_max2, x_C_0spin1, x_C_0spin2, F_A1, M_R1, F_A2, M_R2, M_x_R_ist, M_x_R_delta, F_L, r_delta
    
    #for loop provides sufficient convergence for the equilibrium of forces from normal and tangential forces
    for ii in range(8):
        #Force equilibirum
        F_N1 = 0.5*((-F_ax - np.cos(alpha_rad)*(F_T1+F_T2))/(np.sin(alpha_rad)) + (-F_rad + np.sin(alpha_rad)*(F_T1-F_T2))/(np.cos(alpha_rad)))
        F_N2 = 0.5*((F_ax + np.cos(alpha_rad)*(F_T1+F_T2))/(np.sin(alpha_rad)) + (-F_rad + np.sin(alpha_rad)*(F_T1-F_T2))/(np.cos(alpha_rad)))
        
        #Abort when one contact of the track roller is lifted off the guideway
        if F_N1 < 0 or F_N2 < 0:
            print("--- Caution: Only one contact patch ! ---")
            print("--- Exit Program ---")
            print()
            sys.exit(1)
        
        #Call Hertz.pressung
        a1, b1, p_max1 = Hertz.pressung(F_N1, sumk, E_ers, xi_IST, eta_IST)
        a2, b2, p_max2 = Hertz.pressung(F_N2, sumk, E_ers, xi_IST, eta_IST)
        
        #Calculate the x_C component of the spin center for both contacts
        x_C_0spin1 = r_delta0/np.sin(alpha_rad)
        x_C_0spin2 = -x_C_0spin1
        
        #Call Kinematik.kontaktkraefte mwith short return      
        F_T1, F_A1, M_R1 = Kinematik.kontaktkraefte(n_x_C, a1, b1, p_max1, x_C_0spin1, alpha_rad1, mu, r_KonturFuehrung, r_Abroll, False)
        F_T2, F_A2, M_R2 = Kinematik.kontaktkraefte(n_x_C, a2, b2, p_max2, x_C_0spin2, alpha_rad2, mu, r_KonturFuehrung, r_Abroll, False)
        
        #Toruqe equilibrium in axial direction (x)
        M_x_R_ist = -(M_R1-M_R2)*np.sin(alpha_rad) - (F_A1+F_A2)*r_Abroll
        
        #Force equilibirum in direction of rolling (y)
        F_L = -F_A1-F_A2
        
        #r_delta is a global variable, so that it can be accessed outside of the optimierung_kraefte function as well
        r_delta = r_delta0
        
        #With the function optimierung_kraefte the root of M_x_R_delta is to be found
        M_x_R_delta = M_x_R_ist - M_x_R_soll
        
    return M_x_R_delta


"---------Start of the simulation procedure-------------"

"Load parameters"
#Load all parameters from "Input_and_Run_Sim.py"
n_x_C,v_D,F_Motor,F_ax,F_rad,mu,Material_rolle,Material_Schiene,alpha_deg,A_T,r_KonturRolle,r_KruemmungFuehrung,r_KonturFuehrung,Materialspannungsberechnung=parameter()

#Material der track roller   
#Taking the information about the material from "Materialskript.py"
E_Rolle,nu_Rolle,R_m_Rolle,sigma_D_Rolle,N_Rm_Rolle,N_k_Rolle=Materialskript.Auswahl(Material_rolle)

#Material guideway
#Taking the information about the material from "Materialskript.py"
E_Fuehrung,nu_Fuehrung,R_m_Schiene,sigma_D_Schiene,N_Rm_Schiene,N_k_Schiene=Materialskript.Auswahl(Material_Schiene)


"Calculation"
#Call Hertz.parameter
E_ers, sumk, xi_IST, eta_IST, k_H, n_H, alpha_rad, r_Abroll, h_N, h_T = Hertz.parameter(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, alpha_deg, A_T, r_KonturRolle, r_KruemmungFuehrung, r_KonturFuehrung)

#alpha for both contacts 1 (in rolling direction right) and 2 (in rolling direction left)
alpha_rad1 = alpha_rad
alpha_rad2 = -alpha_rad

#Maximum load at which the roller spins completely
F_L_max_abs = abs(mu*F_rad/np.cos(alpha_rad))


#starting conditions for optimierung_kraefte:
M_x_R_soll = F_Motor*r_Abroll   # Required roller torque in axial direction
F_T1 = 0                        # The tangential forces on both contacts are 
F_T2 = 0                        # set to zero as a starting condition


try:        # If the calculation inside "try:" fails, slip=100% -> it continues at "except ValueError:".
    print("----------------------")
    print("Force optimization")
    print("----------------------")
    #Call optimierung_kraefte
    optimize.brentq(optimierung_kraefte, -r_Abroll*2/3, r_Abroll*2/3, maxiter=250)
    
    print("r_delta = ", r_delta)
    print("M_x_R_delta = ", M_x_R_delta)
    
    #Calculation and output of slip and friction torque
    Schlupf = r_delta/r_Abroll
    Reibmoment = (M_R1-M_R2)*np.sin(alpha_rad)/1000     #Frictional torque with respect to the roller axis in Nm
    print()
    print("Slip = ", round(Schlupf,5))
    print("Frictional Torque M_x_R = ", round(Reibmoment,5), "Nm")
    print()
    
    #Output in the console gives feedback on the calculation
    if M_x_R_soll != 0:
        M_x_R_fehler = abs(M_x_R_delta/M_x_R_soll)*100  #in percent
        if M_x_R_fehler > 5:
            print("--- Caution: Calculation of the contact forces is not accurate! ---")
        else:
            print("Calculation of the contact forces succesfull")
        print("Error: ", round(M_x_R_fehler,3),"%")
    
    print()
    print("----------------------")
    
    
    
    "Plots"
    
    showPlots = True                    # parameter which determines the return type of the function Kinematik.kontaktkraefte
    #Call Kinematik.kontaktkraefte with detailed return
    F_T1, F_A1, x_C_linsp1, y_C_linsp1, x_C_mesh1, y_C_mesh1, v_rel_x_C1, v_rel_y_C1, v_rel_res1, v_rel_x_C_dir1, v_rel_y_C_dir1, v_norm_x_C1, v_norm_y_C1, v_norm_res1, p1, tau_res1, tau_x_C1, tau_y_C1, dP2dxdy1, dPdy1, P1, dW2dxdy1, dWdy1, W1, M_R1, y_C_0spin1 = Kinematik.kontaktkraefte(n_x_C, a1, b1, p_max1, x_C_0spin1, alpha_rad1, mu, r_KonturFuehrung, r_Abroll, showPlots)
    F_T2, F_A2, x_C_linsp2, y_C_linsp2, x_C_mesh2, y_C_mesh2, v_rel_x_C2, v_rel_y_C2, v_rel_res2, v_rel_x_C_dir2, v_rel_y_C_dir2, v_norm_x_C2, v_norm_y_C2, v_norm_res2, p2, tau_res2, tau_x_C2, tau_y_C2, dP2dxdy2, dPdy2, P2, dW2dxdy2, dWdy2, W2, M_R2, y_C_0spin2 = Kinematik.kontaktkraefte(n_x_C, a2, b2, p_max2, x_C_0spin2, alpha_rad2, mu, r_KonturFuehrung, r_Abroll, showPlots)
    
    #Call Plots.multiplot
    Plots.multiplot(F_N1, F_T1, F_A1, x_C_linsp1, y_C_linsp1, x_C_mesh1, y_C_mesh1, v_rel_x_C1, v_rel_y_C1, v_rel_res1, v_rel_x_C_dir1, v_rel_y_C_dir1, p1, tau_res1, tau_x_C1, tau_y_C1, dP2dxdy1, dPdy1, P1, dW2dxdy1, dWdy1, W1, F_N2, F_T2, F_A2, x_C_linsp2, y_C_linsp2, x_C_mesh2, y_C_mesh2, v_rel_x_C2, v_rel_y_C2, v_rel_res2, v_rel_x_C_dir2, v_rel_y_C_dir2, p2, tau_res2, tau_x_C2, tau_y_C2, dP2dxdy2, dPdy2, P2, dW2dxdy2, dWdy2, W2, a1, b1, a2, b2, x_C_0spin1, y_C_0spin1, x_C_0spin2, y_C_0spin2)
    
    #Plotting the forces and contact geometry
    Plots.print_results(F_N1, F_T1, F_A1, F_N2, F_T2, F_A2, h_N, h_T,a1,b1,a2,b2)
   
    
    "Sub-surface-stress calculation & life-time prediction"    
    if Materialspannungsberechnung == "on":
        print()
        print("----------------------")
        print("Start Calculation of max. v.Mises subsurface-stress")
        print("----------------------")
        print("Contact 1")
        #Search global stress maximum: Call Potentialtheorie.vergleichsspannung_max
        vGEH_max, x_opt, y_opt, z_opt = Potentialtheorie.vergleichsspannung_max(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a1, b1, p1, tau_x_C1, tau_y_C1, x_C_linsp1, y_C_linsp1)
        #Calculate the material stresses required for the plots: Call Potentialtheorie.vergleichsspannung_plot
        x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy = Potentialtheorie.vergleichsspannung_plot(vGEH_max, x_opt, y_opt, z_opt, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a1, b1, p1, tau_x_C1, tau_y_C1, x_C_linsp1, y_C_linsp1)
        
        # Calculation of the maximum tolerable number of oscillations for the track roller for contact 1
        N1_max_Rolle,t1_max_Rolle=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m_Rolle,sigma_D_Rolle,N_Rm_Rolle,N_k_Rolle,v_D,r_Abroll)

        # Calculation of the maximum tolerable number of oscillations for the guideway for contact 1
        N1_max_Schiene,t1_max_Schiene=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m_Schiene,sigma_D_Schiene,N_Rm_Schiene,N_k_Schiene,v_D,r_KruemmungFuehrung)         

        #Call Plots.plot_vergleichsspannung
        Plots.plot_vergleichsspannung(x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy, x_opt, y_opt, z_opt)
        print("----------------------")
        print()
        print("Contact 2")
        #Search global stress maximum: Call Potentialtheorie.vergleichsspannung_max
        vGEH_max, x_opt, y_opt, z_opt = Potentialtheorie.vergleichsspannung_max(E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a2, b2, p2, tau_x_C2, tau_y_C2, x_C_linsp2, y_C_linsp2)
        #Calculate the material stresses required for the plots: Call Potentialtheorie.vergleichsspannung_plot
        x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy = Potentialtheorie.vergleichsspannung_plot(vGEH_max, x_opt, y_opt, z_opt, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, a2, b2, p2, tau_x_C2, tau_y_C2, x_C_linsp2, y_C_linsp2)
         
        # Calculation of the maximum tolerable number of oscillations for the track roller for contact 2     
        N2_max_Rolle,t2_max_Rolle=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m_Rolle,sigma_D_Rolle,N_Rm_Rolle,N_k_Rolle,v_D,r_Abroll)

        # Calculation of the maximum tolerable number of oscillations for the guideway for contact 2        
        N2_max_Schiene,t2_max_Schiene=Woehlerlinie.SchwingspielzahlMAX(vGEH_max,R_m_Schiene,sigma_D_Schiene,N_Rm_Schiene,N_k_Schiene,v_D,r_KruemmungFuehrung)              
        
        #Call Plots.plot_vergleichsspannung
        Plots.plot_vergleichsspannung(x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy, x_opt, y_opt, z_opt)
        print("----------------------")

#If slip = 100%,  "try:" fails with"ValueError" -> Continues here:
except ValueError:
    print("--- Caution: Driving force cannot be transmitted ! ---")
    print("--- Instead calculating special case of 100% slip ---")
    #Call Spezialfall.durchdrehen
    Spezialfall.durchdrehen(n_x_C, F_Motor, F_ax, F_rad, E_Rolle, nu_Rolle, E_Fuehrung, nu_Fuehrung, mu, alpha_deg, A_T, r_KonturRolle, r_KruemmungFuehrung, r_KonturFuehrung, Materialspannungsberechnung)



#%%
"-------Post-Processing lifetime-prediction------"
N_Rolle_max_arr=np.array([N1_max_Rolle,N2_max_Rolle])
N_Rolle_max=np.amin(N_Rolle_max_arr)

t_Rolle_max_arr=np.array([t1_max_Rolle,t2_max_Rolle])
t_Rolle_max=np.amin(t_Rolle_max_arr)

print("")
print("----------------------")
print("Prediction")
print("----------------------")
print("Mean velocity = ",v_D," m/min")
print("")

print("First failure track roller")
print("Overrollings = ",round(N_Rolle_max,2))
print("Duration = ",round(t_Rolle_max,2), " h")

N_Schiene_max_arr=np.array([N1_max_Schiene,N2_max_Schiene])
N_Schiene_max=np.amin(N_Schiene_max_arr)

t_Schiene_max_arr=np.array([t1_max_Schiene,t2_max_Schiene])
t_Schiene_max=np.amin(t_Schiene_max_arr)

print("")
print("First failure guideway")
print("Overrollings = ",round(N_Schiene_max,2))
print("Duration = ",round(t_Schiene_max,2), " h")

print("---------END-----------")


