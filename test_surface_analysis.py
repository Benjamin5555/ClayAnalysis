"""
    Provides a simple test of the surface analysis functionality working by reanalysing the 
    surface at each step and checking that the same number of surface atoms are detected at each 
    time step, as we would expect for a clay that is more or less static as in this system. 

"""


import MDAnalysis as mda
import nglview as nv
from nglview.datafiles import PDB, XTC
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from clayAnalysis import ClayAnalysis 
from MDAnalysis import transformations
from MDAnalysis import analysis
from MDAnalysis.analysis import lineardensity


def num_test(lower_surf,upper_surf):
    ag_upper = upper_surf[0]
    n_up = upper_surf[0].n_atoms 
    for i in upper_surf[1:]:
        ag_upper = ag_upper + i
        n_up = n_up+i.n_atoms

    n_low = lower_surf[0].n_atoms
    ag_lower = lower_surf[0]
    for i in lower_surf[1:]:
        ag_lower = ag_lower + i
        n_low = n_low+i.n_atoms
    
    print("waters:" + str((n_low,n_up)))
    #assert (n_low == n_up)
    return ag_upper,ag_lower

def plot_surface_vs_bulk(u,upper_surf,lower_surf,ag_upper,ag_lower):
    clay = u.select_atoms("resname NON*")
    clay = clay - ag_upper 
    clay = clay - ag_lower

    clay_x = np.transpose(clay.positions)[0]
    clay_y = np.transpose(clay.positions)[1]
    clay_z = np.transpose(clay.positions)[2]

    fig = plt.figure()
    ax = Axes3D(fig, zlabel="z")
    ax.set_zlim(50,140)
    for i in range(len(lower_surf)):
        lower_x = np.transpose(lower_surf[i].positions)[0]
        lower_y = np.transpose(lower_surf[i].positions)[1]
        lower_z = np.transpose(lower_surf[i].positions)[2]
        
        upper_x = np.transpose(upper_surf[i].positions)[0]
        upper_y = np.transpose(upper_surf[i].positions)[1]
        upper_z = np.transpose(upper_surf[i].positions)[2]
        
        ax.scatter(upper_x,upper_y,upper_z)
        ax.scatter(lower_x,lower_y,lower_z)
    
    ax.scatter(upper_x,upper_y,upper_z)
    ax.scatter(lower_x,lower_y,lower_z)
    ax.scatter(clay_x,clay_y,clay_z, s=100, lw = 0, color=[1.,1, 1])
    
    clay = u.select_atoms("resname NON*")
    clay = clay - ag_upper 
    clay = clay - ag_lower

    clay_x = np.transpose(clay.positions)[0]
    clay_y = np.transpose(clay.positions)[1]
    clay_z = np.transpose(clay.positions)[2]

    
    fig = plt.figure()
    ax = Axes3D(fig, zlabel="z")
    ax.set_zlim(0,50)
    for i in range(len(lower_surf)):
        lower_x = np.transpose(lower_surf[i].positions)[0]
        lower_y = np.transpose(lower_surf[i].positions)[1]
        lower_z = np.transpose(lower_surf[i].positions)[2]
        
        upper_x = np.transpose(upper_surf[i].positions)[0]
        upper_y = np.transpose(upper_surf[i].positions)[1]
        upper_z = np.transpose(upper_surf[i].positions)[2]
        
        ax.scatter(upper_x,upper_y,upper_z)
        ax.scatter(lower_x,lower_y,lower_z)
    
    ax.scatter(upper_x,upper_y,upper_z)
    ax.scatter(lower_x,lower_y,lower_z)
    ax.scatter(clay_x,clay_y,clay_z, s=100, lw = 0, color=[1.,1, 1])
    
    plt.show()

#u = mda.Universe("TestFiles/eq_d_spacing_06.tpr","TestFiles/eq_d_spacing_06.trr")
#u = mda.Universe("TestFiles/sim_01.tpr","TestFiles/sim_01.trr")
u = mda.Universe("TestFiles/sim_02/sim_02.tpr","TestFiles/sim_02/sim_02.trr")


if __name__ == "__main__": 
    
    cal = ClayAnalysis(u)
    
    
    
    for du in u.trajectory:
        surfaces = cal.generate_surface_group("waters")
        lower_surf = surfaces[0]
        upper_surf = surfaces[1]
        ag_upper,ag_lower = num_test(lower_surf,upper_surf) 
        #plot_surface_vs_bulk(u,upper_surf,lower_surf,ag_upper,ag_lower)
        assert cal.combine_atomgroups(lower_surf)==ag_lower
    
    
        surfaces = cal.generate_surface_group("mineral")
        lower_surf = surfaces[0]
        upper_surf = surfaces[1]
        #
        ag_upper,ag_lower = num_test(lower_surf,upper_surf)
        #`plot_surface_vs_bulk(u,upper_surf,lower_surf,ag_upper,ag_lower)
    
        
    
