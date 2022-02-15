
"""

Technical Universitiy Darmstadt, Germany
Institute for Product Development and Machine Elements

Defines input parameters and starts the simulation

The simulation will start when running this file.

To have access to the full ouput, you can start the "_main.py" script directly.

"""

def parameter():

    "Forces"
    #The description of the coordinate systems and direction of the forces can be found in the documentation
    #F_Motor > 0 -> Braking 
    #F_Motor < 0 -> Driving
    F_Motor = 30            # Driving Force [N]
    F_ax = 30               # Axial Force [N]        #positive axial load in x-direction 
    F_rad = -200             # Radial Force [N]      #negative radial load in z-direction is equivalent to a pre-load
    
    "Material"
    #Material of the track roller   
    Material_track_roller="1.7225-unvergütet"
    #Material of the guideway
    Material_guideway="1.7225-unvergütet"
    # Coefficient of Friction [-]
    mu = 0.3                
    
    "Geometry"    
    #Geometry of track roller and guideway
    alpha_deg = 30                  # contact angle [°]
    
    A_T = 21.75                     # [mm]
    r_C_TR = -7              # [mm]
    
    r_Curve = 39.5      # radius of the curve of the guideway [mm], for straight guideways -> np.inf
    r_C_G = 6            # [mm]
    
    "Velocity"
    #Mean velocity [m/min] (For life-time prediction)
    v_D = 60 
    
    "Simulation Settings"
    #Number of discrete points in x-direction
    n_x_C = 40     #Values between 40 (min. 20) and 80 (max. 120) are recommanded 
    #Sub-surface stress calculation option
    subsurface_stress = "on"  # If not "on" -> no sub-surface stress calculation

    return n_x_C,v_D,F_Motor,F_ax,F_rad,mu,Material_track_roller,Material_guideway,alpha_deg,A_T,r_C_TR,r_Curve,r_C_G,subsurface_stress
            
#Runs the simulation
def run():
    import _main
    _main

run()