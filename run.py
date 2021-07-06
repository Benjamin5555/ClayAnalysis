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
BINS= 500

def shift_clay_to_origin(ts):
    """
    Translate all coordinates to have clay centered.

    Might be better to pass clay atom selection? simpler version for now
    """
    clay = u.select_atoms("resname NON*")
    
    
    ts.positions = ts.positions - clay.center_of_mass()
    return ts

u = mda.Universe("TestFiles/sim_02.tpr","TestFiles/sim_02.trr")

#u = mda.Universe("TestFiles/eq_d_spacing_06.tpr","TestFiles/eq_d_spacing_06.trr")
cal = ClayAnalysis(u)
clay = u.select_atoms("resname NON*")

u.trajectory.add_transformations(*[shift_clay_to_origin,mda.transformations.wrap(u.atoms)])

#w = nv.show_mdanalysis(u)
#w


####################Clay density mostly from DynDol############
#box_dims = cal.get_box_dim()
##print(np.unique(u.atoms.names))
#
#partial_densities = []
#print("Determining partial densities:")
#for sel in ["O[!W]*","ST*","AT*"]:#np.unique(u.atoms.names): 
#    sel= "name "+sel
#    partial_densities.append(cal.get_partial_density(box_dims, sel, bins=100))
#
#
#partial_densities.append(cal.get_partial_density(box_dims, "resname SOL", bins=100))
#
#
#print("\n\nClay density time average along z\n\n ")
#cal.plot_average_density_t(partial_densities)
#
#plt.legend()
#
#
########################
#
#"""
#Clay surface determination
#
#(assumes surface mostly still over time, but could change bounds of surface every time step
#"""

for du in u.trajectory:
    surfaces = cal.generate_surface_group("mineral")
    lower_surf = surfaces[0]
    upper_surf = surfaces[1]
    
#    fig = plt.figure()
#    ax = Axes3D(fig, zlabel="z")
#    for i in range(len(lower_surf)):
#        lower_x = np.transpose(lower_surf[i].positions)[0]
#        lower_y = np.transpose(lower_surf[i].positions)[1]
#        lower_z = np.transpose(lower_surf[i].positions)[2]
#        
#        upper_x = np.transpose(upper_surf[i].positions)[0]
#        upper_y = np.transpose(upper_surf[i].positions)[1]
#        upper_z = np.transpose(upper_surf[i].positions)[2]
#        
#        ax.scatter(upper_x,upper_y,upper_z)
#        ax.scatter(lower_x,lower_y,lower_z)
    #plt.show()    

    ag_upper = upper_surf[0]
    n_up = upper_surf[0].n_atoms 
    for i in upper_surf:
        ag_upper = ag_upper + i
        n_up = n_up+i.n_atoms

    n_low = lower_surf[0].n_atoms
    ag_lower = lower_surf[0]
    for i in lower_surf:
        ag_lower = ag_lower + i
        n_low = n_low+i.n_atoms
    
    clay = clay - ag_upper 
    clay = clay - ag_lower
    print(n_low,n_up)
    #assert (n_low == n_up)
    print("\n\nComparison of surface to whole clay")
    clay_x = np.transpose(clay.positions)[0]
    clay_y = np.transpose(clay.positions)[1]
    clay_z = np.transpose(clay.positions)[2]
   
    
    fig = plt.figure()
    ax = Axes3D(fig, zlabel="z")
    ax.set_zlim(110,113)
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
    
    
    
    fig = plt.figure()
    ax = Axes3D(fig, zlabel="z")
    ax.set_zlim(20,24)
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
    






















###########################################################
#
#print("\n\nadsorption of Mg to Si by looking at distance between (only looking at 'top' surface ")
#r_c = 5 # 5 angstrom cutoff distance
#
#layer = u.select_atoms("prop z >= 87.25 and prop z < 89.36") #Surface determination doesn't actually
#                                                             #work with this yet in differnt forms
#t_surface_Si = layer.atoms.select_atoms("name S*")
#test_ads_names =  [u.select_atoms("resname Mg").names[0]]
#t_surface_atoms_ids = t_surface_Si.indices
#
#adsorbed = cal.find_adsorbed(t_surface_atoms_ids,test_ads_names,r_c)
#
##print(adsorbed)
#a_x = np.transpose(adsorbed.positions)[0]
#a_y = np.transpose(adsorbed.positions)[1]
#a_z = np.transpose(adsorbed.positions)[2]
#
#t_surface_Si_x = np.transpose(t_surface_Si.positions)[0]
#t_surface_Si_y = np.transpose(t_surface_Si.positions)[1]
#t_surface_Si_z = np.transpose(t_surface_Si.positions)[2]
#
#layer_x = np.transpose(layer.positions)[0]
#layer_y = np.transpose(layer.positions)[1]
#layer_z = np.transpose(layer.positions)[2]
#
#fig = plt.figure()
#ax = Axes3D(fig, zlabel="z")
#
#ax.scatter(a_x,a_y,a_z)
#
#ax.scatter(t_surface_Si_x,t_surface_Si_y,t_surface_Si_z)
#
#plt.show()

minima = get_minima_coords(density[-1][1],density[-1][0])

#Want to find the z axis boundary of the lower clay basal surface
lower_surface_bounds = minima[:2]
#Want to find the z axis boundary of the other clay basal surface
upper_surface_bounds = minima[-2:]

