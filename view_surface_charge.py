"""
    Provides a very basic visulisation of charge on the surface of the clay via a 2D histogram
    
    Probably a better idea to just use VMD

"""

from clayAnalysis import ClayAnalysis 
import numpy as np
import MDAnalysis as mda
if __name__ == "__main__":
    u = mda.Universe("TestFiles/control.tpr","TestFiles/control.trr")
    
    cal = ClayAnalysis(u)
    
    surfaces = cal.generate_surface_group("waters")
    upper_surf = surfaces[0]
    lower_surf = surfaces[1]
    
    ag_lower = lower_surf[0]
    
    
    
    
    ag_upper = upper_surf[0]
    for i in upper_surf[1:]:
        ag_upper = ag_upper + i
    
    surfaces = cal.generate_surface_group("mineral")
    upper_surf = surfaces[0]
    for i in upper_surf[1:]:
        ag_upper = ag_upper + i
    
    


#cal.plot_surface_charge(ag_upper,2)

