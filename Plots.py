"""

Technical Universitiy Darmstadt, Germany
Institute for Product Development and Machine Elements

Dynamic Contact Simulation for Locating Bearing Track Rollers

Plots and Results

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle


def multiplot(F_N1, F_T1, F_A1, x_C_linsp1, y_C_linsp1, x_C_mesh1, y_C_mesh1, v_rel_x_C1, v_rel_y_C1, v_rel_res1, v_rel_x_C_dir1, v_rel_y_C_dir1, p1, tau_res1, tau_x_C1, tau_y_C1, dP2dxdy1, dPdy1, P1, dW2dxdy1, dWdy1, W1, F_N2, F_T2, F_A2, x_C_linsp2, y_C_linsp2, x_C_mesh2, y_C_mesh2, v_rel_x_C2, v_rel_y_C2, v_rel_res2, v_rel_x_C_dir2, v_rel_y_C_dir2, p2, tau_res2, tau_x_C2, tau_y_C2, dP2dxdy2, dPdy2, P2, dW2dxdy2, dWdy2, W2, a1, b1, a2, b2, x_C_0spin1, y_C_0spin1, x_C_0spin2, y_C_0spin2):
    """Creates a multiplot visualizing the tangential stresses and the frictional power distributions"""
    
    fig, axs = plt.subplots(2,2, figsize = (18,12), sharex = True, sharey = "row", gridspec_kw={'hspace': 0.25, 'wspace': 0.1})
    
    a_max = np.maximum(a1,a2)
    b_max = np.maximum(b1,b2)
    
    q2 = axs[0,0].quiver(x_C_mesh2, y_C_mesh2, tau_x_C2, tau_y_C2, tau_res2, units="xy", pivot="tip",
                         cmap="viridis"
                         )
    axs[0,0].quiverkey(q2, 0.45, 0.9, 200, r'$200 MPa$', labelpos='E',
                       coordinates='figure')
    # axs[0,0].scatter(x_C_mesh2, y_C_mesh2, color='0.2', s=1)
    axs[0,0].set_xlim(-a_max*1.2, a_max*1.2)
    axs[0,0].set_ylim(-b_max*1.1, b_max*1.1)
    #axs[0,0].set_xlabel('$x_{C2}$ in mm',fontsize=20)
    #axs[0,0].set_ylabel('$y_{C2}$ in mm',fontsize=20)
    axs[0,0].set_ylabel('$y_{C}$ in mm',fontsize=20)
    axs[0,0].set_title("Shear stress", loc = "left",size=16)
    for ii in np.linspace(0,1,201):
        ellipse = Ellipse((0,0), 2*a2 + ii*a2, 2*b2 + ii*b2,edgecolor = "w", linewidth = 2, fill = False)
        axs[0,0].add_patch(ellipse)
    ellipse = Ellipse((0,0), 2*a2, 2*b2,edgecolor = "black", linewidth = 2, fill = False)
    axs[0,0].add_patch(ellipse)
    circle = Circle((x_C_0spin2, y_C_0spin2), radius = a_max/40, color = "black")
    axs[0,0].add_patch(circle)
    axs[0,0].tick_params(axis='y', labelsize=16)
    
    q1 = axs[0,1].quiver(x_C_mesh1, y_C_mesh1, tau_x_C1, tau_y_C1, tau_res1, units="xy", pivot="tip")
    axs[0,1].quiverkey(q1, 0.85, 0.9, 200, r'$200 MPa$', labelpos='E',
                       coordinates='figure')
    # axs[0,1].scatter(x_C_mesh1, y_C_mesh1, color='0.2', s=1)
    #axs[0,1].set_xlabel('$x_{C1}$ in mm')
    #axs[0,1].set_ylabel('$y_{C1}$ in mm')
    axs[0,1].set_title("Shear stress", loc = "left",size=16)
    for ii in np.linspace(0,1,201):
        ellipse = Ellipse((0,0), 2*a1 + ii*a1, 2*b1 + ii*b1,edgecolor = "w", linewidth = 2, fill = False)
        axs[0,1].add_patch(ellipse)
    ellipse = Ellipse((0,0), 2*a1, 2*b1,edgecolor = "black", linewidth = 2, fill = False)
    axs[0,1].add_patch(ellipse)
    circle = Circle((x_C_0spin1, y_C_0spin1), radius = a_max/40, color = "black")
    axs[0,1].add_patch(circle)
    
    
    #plt.xticks(fontsize=16)
    axs[1,0].plot(x_C_linsp2, dWdy2,color="black",linewidth=3)
    axs[1,0].set_xlabel('$x_{C2}$ in mm',fontsize=20)
    axs[1,0].set_ylabel(r'Rel. frictional energy distr. in $\frac{kJ}{m^2}$',fontsize=20)
    axs[1,0].set_title("Fricional energy per rotation: " + str(round(W2, 3)) + r"$\frac{J}{m}$", loc = "left",size=16)
    #axs[1,0].xticks(fontsize=16)
    axs[1,0].grid()
    axs[1,0].tick_params(axis='x', labelsize=16)
    axs[1,0].tick_params(axis='y', labelsize=16)
    
    #plt.xticks(fontsize=16)
    axs[1,1].plot(x_C_linsp1, dWdy1,color="black",linewidth=3)
    axs[1,1].set_xlabel('$x_{C1}$ in mm',fontsize=20)
    #axs[1,1].set_ylabel(r'Relative Reibenergieverteilung in $\frac{kJ}{m^2}$')
    axs[1,1].set_title("Fricional energy per rotation: " + str(round(W1, 3)) + r"$\frac{J}{m}$", loc = "left",size=16)
    axs[1,1].grid()
    axs[1,1].tick_params(axis='x', labelsize=16)
    

    
def print_results(F_N1, F_T1, F_A1, F_N2, F_T2, F_A2, h_N, h_T,a1,b1,a2,b2):
    """Outputs the forces, moments and the sizes of the ellipse in the console"""
    
    print()
    print("Calculated contact forces:")
    print("Normal forces:\t\tF_N2 =", round(F_N2,2), "N\t\tF_N1 =", round(F_N1,2), "N")
    print("Tangential forces:\tF_T2 =", round(F_T2,2), "N\t\tF_T1 =", round(F_T1,2), "N")
    print("Driving forces:\t\tF_A2 =", round(F_A2,2), "N\t\tF_A1 =", round(F_A1,2), "N")
    print()
    M_y_R = h_N*(F_N1-F_N2) + h_T*(F_T1+F_T2)
    print("Tilting torque M_y_R =", round(M_y_R/1000,2), "Nm")
    print()
    print("Calculated contact area:")
    print("Semi-axis a:\t\ta2 =", round(a2,2), "mm\t\ta1 =", round(a1,2), "mm")
    print("Semi-axis b:\t\ta2 =", round(b2,2), "mm\t\ta1 =", round(b1,2), "mm")


def plot_vergleichsspannung(x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy, x_opt, y_opt, z_opt):
    """Creates plots showing the stress distributions inside the material in two sectional planes"""
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (18,4),sharey = True)
    
    r_marker = abs(x_linsp[0])/100
    
    con1 = ax1.contour(x_linsp, z_linsp, vGEH_zx, 12)
    ax1.set_xlabel("x in mm",fontsize=20)
    ax1.set_ylabel("z in mm",fontsize=20)
    ax1.set_title("DEH in MPa" % (round(x_opt,4)), loc = "left",size=16)
    ax1.set_title("$y_{vGEH_{max}}$ = %gmm" % (round(y_opt,4)), loc = "right",size=16)
    ax1.clabel(con1, inline=1, fontsize=10, fmt='%1.0f')
    #fig.colorbar(con1, ax=ax1)x
    max_marker1 = Circle((x_opt, z_opt), radius = r_marker, color = "k")
    ax1.add_patch(max_marker1)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)

    con2 = ax2.contour(y_linsp, z_linsp, vGEH_zy, con1.levels)
    ax2.set_xlabel("y in mm",fontsize=20)
    #ax2.set_ylabel("z-Position in mm")
    ax2.set_title("DEH in MPa" % (round(x_opt,4)), loc = "left",size=16)
    ax2.set_title("$x_{vGEH_{max}}$ = %gmm" % (round(x_opt,4)), loc = "right",size=16)
    ax2.clabel(con2, inline=1, fontsize=10, fmt='%1.0f')
    #fig.colorbar(con2, ax=ax2)
    max_marker2 = Circle((y_opt, z_opt), radius = r_marker, color = "k")
    ax2.add_patch(max_marker2)
    ax2.tick_params(axis='x', labelsize=16)

