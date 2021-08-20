from clayAnalysis import ClayAnalysis 
import MDAnalysis as mda
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

if __name__ == "__main__": 
    #Create universe object from file
    u = mda.Universe("TestFiles/control.tpr","TestFiles/control.trr")
    
    #Create clay analysis object, linking to the universe under consideration
    cal = ClayAnalysis(u)
    
    
    #Caluculate basal surfaces
    clay_mineral_surf = cal.generate_surface_group("mineral")
    
    #Returned list has list of atomgroups in at lower basal surface (transformed so actually top 
    #surface) at 0 index and list of atomgroups at upper surface in index 1
    
    clay_mineral_lower_surf = clay_mineral_surf[0]
    clay_mineral_upper_surf = clay_mineral_surf[1]
    
    #Surface 'minerals' ST & AT and surface waters separatly found 
    clay_waters_surf = cal.generate_surface_group("waters")
    clay_waters_lower_surf = clay_waters_surf[0]
    clay_waters_upper_surf = clay_waters_surf[1]
    
    
    
    
    #Can then work with these groups, combine etc e.g. plotting surface:
    fig = plt.figure()
    ax = Axes3D(fig, zlabel="z")
    for i in range(len(clay_mineral_lower_surf)):
        lower_x = np.transpose(clay_mineral_lower_surf[i].positions)[0]
        lower_y = np.transpose(clay_mineral_lower_surf[i].positions)[1]
        lower_z = np.transpose(clay_mineral_lower_surf[i].positions)[2]
        
        upper_x = np.transpose(clay_mineral_upper_surf[i].positions)[0]
        upper_y = np.transpose(clay_mineral_upper_surf[i].positions)[1]
        upper_z = np.transpose(clay_mineral_upper_surf[i].positions)[2]
     
        ax.scatter(upper_x,upper_y,upper_z)
        ax.scatter(lower_x,lower_y,lower_z)
    
    for i in range(len(clay_waters_lower_surf)):
    
        lower_x = np.transpose(clay_waters_lower_surf[i].positions)[0]
        lower_y = np.transpose(clay_waters_lower_surf[i].positions)[1]
        lower_z = np.transpose(clay_waters_lower_surf[i].positions)[2]
        upper_x = np.transpose(clay_waters_upper_surf[i].positions)[0]
        upper_y = np.transpose(clay_waters_upper_surf[i].positions)[1]
        upper_z = np.transpose(clay_waters_upper_surf[i].positions)[2]
       
        ax.scatter(upper_x,upper_y,upper_z)
        ax.scatter(lower_x,lower_y,lower_z)
    
    ax.legend()
    plt.show()    
    
