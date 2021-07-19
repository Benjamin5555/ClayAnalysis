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

def plot_group(grp):
    
    fig = plt.figure()
    ax = Axes3D(fig, zlabel="z")

    upper_x = np.transpose(grp.positions)[0]
    upper_y = np.transpose(grp.positions)[1]
    upper_z = np.transpose(grp.positions)[2]


    ax.scatter(upper_x,upper_y,upper_z)

    plt.show()






##################TEST 2
print("Test2")
u = mda.Universe("TestFiles/test_sys.pdb","TestFiles/test_sys.pdb")

#u = mda.Universe("TestFiles/sim_02/sim_02.tpr","TestFiles/test_sys.pdb")

#plot_group(u.atoms)
cal = ClayAnalysis(u)


layer = u.atoms.select_atoms("not resname Ca") #and (name AT* or name ST*)
ads = u.atoms.select_atoms("resname Ca") #and (name AT* or name ST*)")

surface_atoms_ids = layer.atoms.indices
test_ads_names = ads.atoms.resnames
resnames = []
for sid in surface_atoms_ids:
    resnames.append(u.atoms.select_atoms("index "+str(sid)).resnames)

times, stats= cal.adsorption_times(surface_atoms_ids,test_ads_names,1)

print(stats)
print(times)

#All desorbed
assert stats[0] == [0,0,0,0,0]

#Single Adsorption
assert stats[1] ==[1,0,0,1,1]

#Single desorbs 
assert stats[2] == [0,1,0,0,1]

#All adsorb
assert stats[3] == [17,0,0,17,18]

#All Desorb bar one and adsorb to same ion
assert stats[4] == [16,16,1,17,18+16]

#All desorb
assert stats[5] == [0,17,0,0,18+16]

#All adsorb to two surface ions at same time
assert stats[6] == [34,0,0,34,18+16+34]

#All desorb
assert stats[7] == [0,34,0,0,18+16+34]

assert len(np.where(np.array(times)==1)[0])==18+16+34-1
assert len(np.where(np.array(times)==2)[0])==1


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
