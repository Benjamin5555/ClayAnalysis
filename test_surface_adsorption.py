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
assert stats[2] ==[1,0,1,2,2]

#Single desorbs 
assert stats[3] == [0,1,1,1,2]
assert stats[4] == [0,1,0,0,2]

#All adsorb
assert stats[5] == [17,0,0,17,19]

#All Desorb bar one and adsorb to same ion
assert stats[6] == [16,16,1,17,19+16]


#All desorb
assert stats[7] == [0,17,0,0,19+16]

#All adsorb to two surface ions at same time
assert stats[8] == [34,0,0,34,19+16+34]

#All desorb
assert stats[9] == [0, 34, 0, 0,19+16+34]

#ALL adsorb to one again 
assert stats[10] == [17, 0, 0,17,19+16+34+17]

#8 desorb from one 
assert stats[11] == [ 0, 8, 17-8,17-8,86]

#4 more desorb
assert stats[12] == [0,5,17-8-5,17-8-5,86]

#ALL BUT 1 desorb
assert stats[13] == [ 0, 3, 1, 1,86]

#CHECK TIMES OUTPUT, expect mostly 1 timestep adsorptions but 10 for 2 steps 



assert len(np.where(np.array(times)==1)[0])==75
assert len(np.where(np.array(times)==2)[0])==6 #One ion remains joined which is not added to the times

assert len(np.where(np.array(times)==3)[0])==4



