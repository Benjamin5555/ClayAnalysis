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



u = mda.Universe("TestFiles/eq_d_spacing_06.tpr","TestFiles/eq_d_spacing_06.trr")
#u = mda.Universe("TestFiles/sim_02.tpr","TestFiles/sim_02.trr")

cal = ClayAnalysis(u)

layer = u.select_atoms("prop z >= 87.25 and prop z < 89.36") #Surface determination doesn't actually
                                                             #work with this yet in differnt forms

r_c = 5 # 5 angstrom cutoff distance
t_surface_Si = layer.atoms.select_atoms("name S*")
test_ads_names =  [u.select_atoms("resname Mg").names[0]]
t_surface_atoms_ids = t_surface_Si.indices


adsorbed = cal.find_adsorbed_as_agroup(t_surface_atoms_ids,test_ads_names,r_c)
print(adsorbed)
#fig = plt.figure()
#ax = Axes3D(fig, zlabel="z")
#t_surface_Si_x = np.transpose(t_surface_Si.positions)[0]
#t_surface_Si_y = np.transpose(t_surface_Si.positions)[1]
#t_surface_Si_z = np.transpose(t_surface_Si.positions)[2]
#
#for ads in adsorbed.keys:
#    a_x = np.transpose(adsorbed[ads].positions)[0]
#    a_y = np.transpose(adsorbed[ads].positions)[1]
#    a_z = np.transpose(adsorbed[ads].positions)[2]
#        
#    
#    
#    ax.scatter(a_x,a_y,a_z)
#    
#ax.scatter(t_surface_Si_x,t_surface_Si_y,t_surface_Si_z)
#
#plt.show()
