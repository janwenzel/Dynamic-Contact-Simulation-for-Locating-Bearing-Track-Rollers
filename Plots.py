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
    """Erstellt einen Multiplot, der die Tangentialspannungs- und die Leistungsverteilungen visualisiert"""
    
    fig, axs = plt.subplots(2,2, figsize = (12,8), sharex = True, sharey = "row", gridspec_kw={'hspace': 0.25, 'wspace': 0.1})
    
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
    axs[0,0].set_xlabel('$x_{C2}$ in mm')
    axs[0,0].set_ylabel('$y_{C2}$ in mm')
    axs[0,0].set_title("Tangentialspannung", loc = "left")
    for ii in np.linspace(0,1,201):
        ellipse = Ellipse((0,0), 2*a2 + ii*a2, 2*b2 + ii*b2,edgecolor = "w", linewidth = 2, fill = False)
        axs[0,0].add_patch(ellipse)
    ellipse = Ellipse((0,0), 2*a2, 2*b2,edgecolor = "black", linewidth = 2, fill = False)
    axs[0,0].add_patch(ellipse)
    circle = Circle((x_C_0spin2, y_C_0spin2), radius = a_max/40, color = "black")
    axs[0,0].add_patch(circle)
    
    
    q1 = axs[0,1].quiver(x_C_mesh1, y_C_mesh1, tau_x_C1, tau_y_C1, tau_res1, units="xy", pivot="tip")
    axs[0,1].quiverkey(q1, 0.87, 0.9, 200, r'$200 MPa$', labelpos='E',
                       coordinates='figure')
    # axs[0,1].scatter(x_C_mesh1, y_C_mesh1, color='0.2', s=1)
    axs[0,1].set_xlabel('$x_{C1}$ in mm')
    axs[0,1].set_ylabel('$y_{C1}$ in mm')
    axs[0,1].set_title("Tangentialspannung", loc = "left")
    for ii in np.linspace(0,1,201):
        ellipse = Ellipse((0,0), 2*a1 + ii*a1, 2*b1 + ii*b1,edgecolor = "w", linewidth = 2, fill = False)
        axs[0,1].add_patch(ellipse)
    ellipse = Ellipse((0,0), 2*a1, 2*b1,edgecolor = "black", linewidth = 2, fill = False)
    axs[0,1].add_patch(ellipse)
    circle = Circle((x_C_0spin1, y_C_0spin1), radius = a_max/40, color = "black")
    axs[0,1].add_patch(circle)
    
    
    
    axs[1,0].plot(x_C_linsp2, dWdy2,color="black",linewidth=3)
    axs[1,0].set_xlabel('$x_{C2}$ in mm')
    axs[1,0].set_ylabel(r'Relative Reibenergieverteilung in $\frac{kJ}{m^2}$')
    axs[1,0].set_title("Reibenergie pro Rollweg: " + str(round(W2, 3)) + r"$\frac{J}{m}$", loc = "left")
    axs[1,0].grid()
    
    
    axs[1,1].plot(x_C_linsp1, dWdy1,color="black",linewidth=3)
    axs[1,1].set_xlabel('$x_{C1}$ in mm')
    axs[1,1].set_ylabel(r'Relative Reibenergieverteilung in $\frac{kJ}{m^2}$')
    axs[1,1].set_title("Reibenergie pro Rollweg: " + str(round(W1, 3)) + r"$\frac{J}{m}$", loc = "left")
    axs[1,1].grid()
    
    
def print_results(F_N1, F_T1, F_A1, F_N2, F_T2, F_A2, h_N, h_T,a1,b1,a2,b2):
    """Gibt die Kräfte, Momenten und Größen der Ellipse in der Konsole aus"""
    
    print()
    print("Berechnete Kontaktkräfte:")
    print("Normalkräfte:\t\tF_N2 =", round(F_N2,2), "N\t\tF_N1 =", round(F_N1,2), "N")
    print("Tangentialkräfte:\tF_T2 =", round(F_T2,2), "N\t\tF_T1 =", round(F_T1,2), "N")
    print("Antriebskräfte:\t\tF_A2 =", round(F_A2,2), "N\t\tF_A1 =", round(F_A1,2), "N")
    print()
    M_y_R = h_N*(F_N1-F_N2) + h_T*(F_T1+F_T2)
    print("Kippmoment M_y_R =", round(M_y_R/1000,2), "Nm")
    print()
    print("Berechnete Kontaktfläche:")
    print("Halbachse a:\t\ta2 =", round(a2,2), "mm\t\ta1 =", round(a1,2), "mm")
    print("Halbachse b:\t\ta2 =", round(b2,2), "mm\t\ta1 =", round(b1,2), "mm")


def plot_vergleichsspannung(x_linsp, y_linsp, z_linsp, vGEH_zx, vGEH_zy, x_opt, y_opt, z_opt):
    """Erstellt Plots, welche die Spannungsverteilungen im Materialinneren in zwei Schnittebenen zeigen"""
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (18,4))
    
    r_marker = abs(x_linsp[0])/100
    
    con1 = ax1.contour(x_linsp, z_linsp, vGEH_zx, 12)
    ax1.set_xlabel("x-Position in mm")
    ax1.set_ylabel("z-Position in mm")
    ax1.set_title("Vergleichsspannungen in MPa" % (round(x_opt,4)), loc = "left")
    ax1.set_title("$y_{vGEH_{max}}$ = %gmm" % (round(y_opt,4)), loc = "right")
    ax1.clabel(con1, inline=1, fontsize=10, fmt='%1.0f')
    #fig.colorbar(con1, ax=ax1)x
    max_marker1 = Circle((x_opt, z_opt), radius = r_marker, color = "k")
    ax1.add_patch(max_marker1)

    con2 = ax2.contour(y_linsp, z_linsp, vGEH_zy, con1.levels)
    ax2.set_xlabel("y-Position in mm")
    ax2.set_ylabel("z-Position in mm")
    ax2.set_title("Vergleichsspannungen in MPa" % (round(x_opt,4)), loc = "left")
    ax2.set_title("$x_{vGEH_{max}}$ = %gmm" % (round(x_opt,4)), loc = "right")
    ax2.clabel(con2, inline=1, fontsize=10, fmt='%1.0f')
    #fig.colorbar(con2, ax=ax2)
    max_marker2 = Circle((y_opt, z_opt), radius = r_marker, color = "k")
    ax2.add_patch(max_marker2)


