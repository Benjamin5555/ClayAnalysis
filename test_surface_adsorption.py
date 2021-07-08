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
import pytest


def plot_group_2d(grp):
    
    upper_x = np.transpose(grp.positions)[0]
    
    upper_z = np.transpose(grp.positions)[2]
        # non_surf = clay.select_atoms("resname "+sel).difference(ag_upper)
        # print(non_surf)

    plt.scatter(upper_x,upper_z)

def plot_group(grp,fig,ax):

    upper_x = np.transpose(grp.positions)[0]
    upper_y = np.transpose(grp.positions)[1]
    upper_z = np.transpose(grp.positions)[2]

        # non_surf = clay.select_atoms("resname "+sel).difference(ag_upper)
        # print(non_surf)

    ax.scatter(upper_x,upper_y,upper_z)




#u = mda.Universe("TestFiles/test_sys.pdb","TestFiles/test_sys.pdb")
print("Bulk Charcteristics check")

u = mda.Universe("TestFiles/test_sys.pdb","TestFiles/test_sys.pdb")

cal = ClayAnalysis(u)

layer = u.select_atoms("resname NON*") 

upper_x = np.transpose(layer.positions)[0]
    
upper_z = np.transpose(layer.positions)[2]
        # non_surf = clay.select_atoms("resname "+sel).difference(ag_upper)
        # print(non_surf)
plt.scatter(upper_x,upper_z)

test_ads_names =  [u.select_atoms("resname Mg").names[0]]
surface_atoms_ids = layer.indices

#fig = plt.figure()
#ax = Axes3D(fig, zlabel="z")
Mgs = u.select_atoms("type Mg")
plot_group_2d(layer)
plot_group_2d(Mgs)
#`plt.show()


times, stats= cal.adsorption_times(surface_atoms_ids,test_ads_names,0,8)

#NO CHANGE/ADSORPTION
assert stats[0] == [0,0,0,0]

#ADSORPTION ONLY
assert stats[1][0] != 0 and stats[1][1] == 0 
#NO CHANGE
assert stats[2][0] == 0 and stats[2][1] == 0 
#ALL DESORB
assert stats[3][2] == stats[2][1] and stats[3][3] == 0
#ALL ADSORB
assert stats[4][0] != 0 and stats[4][1] == 0 
#HALF DESORB
assert stats[5][1] == pytest.approx(stats[5][3],1)
#HALF DESORB, HALF ADSORB
assert stats[6][0] == pytest.approx(stats[6][1],1) 


#adsorbed = cal.find_adsorbed(surface_atoms_ids,test_ads_names,r_c)

#for ts in u.trajectory[1:20]:
#    #print(Mgs.positions)
#    #Mgs.positions[0] = adsorbed[list(adsorbed.keys())[0]].positions[0]
#
#    #Mgs.positions[np.random.randint(Mgs.n_atoms)][2] = adsorbed[list(adsorbed.keys())[0]].positions[0][2]
#    #print("-----------------------")
#    #print(Mgs.positions)
#    plot_group_2d(layer)
#    adsorbed = cal.find_adsorbed(surface_atoms_ids,test_ads_names,r_c)
#    for surf in adsorbed.keys():
#        print(adsorbed[surf].positions)
#        plot_group_2d(adsorbed[surf])
#    plt.show()



#print(adsorbed)
#
#plot_group_2d(layer)
#plot_group_2d(Mgs)
#for surf in adsorbed.keys():
#   plot_group_2d(adsorbed[surf])
#
#
#plt.show()
