from clayAnalysis import ClayAnalysis 
import MDAnalysis as mda
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import nglview as nv
import copy



def plot_group(grp,fig=None,ax=None,label=None):
    if fig ==None:
        fig = plt.figure()
        ax = Axes3D(fig, zlabel="z")

    upper_x = np.transpose(grp.positions)[0]
    upper_y = np.transpose(grp.positions)[1]
    upper_z = np.transpose(grp.positions)[2]

        # non_surf = clay.select_atoms("resname "+sel).difference(ag_upper)
        # print(non_surf)

    ax.scatter(upper_x,upper_y,upper_z)

def move_away_z(ts):
    """
    Translate all coordinates by 2 angstroms up along the Z dimension.
    """
    ts.positions = ts.positions + np.array([0, 0, 100], dtype=np.float32)
    return ts

def set_z_away(na):
    pos_ad = na.atoms.positions.T 
    pos_ad[2] = pos_ad[2]*0
    return pos_ad.T

#Create universe object from file
u = mda.Universe("TestFiles/sim_02/sim_02.tpr","TestFiles/sim_02/sim_02.trr")


#Create clay analysis object, linking to the universe under consideration
cal = ClayAnalysis(u)


#Caluculate basal surfaces
clay_mineral_surf = cal.generate_surface_group("mineral")

#Returned list has list of atomgroups in at lower basal surface (transformed so actually top 
#surface) at 0 index and list of atomgroups at upper surface in index 1

cmls = cal.combine_atomgroups(clay_mineral_surf[0])

#Surface 'minerals' ST & AT and surface waters separatly found 
clay_waters_surf = cal.generate_surface_group("waters")
cwls =  cal.combine_atomgroups(clay_waters_surf[0])

cw_orig_pos =np.array( cwls.atoms.positions)
#cwls_away = cwls.add_transformation(move_away_z)
print("generating test data")
r_c = 2  
with mda.Writer("test_comp.pdb",u.atoms.n_atoms) as W:
    
    #NOTHING ADSORBED FIRST FRAME
    frame1 = cmls
    na = u.atoms - frame1
    na.atoms.positions = set_z_away(na)
    print("1:")
    print(frame1.atoms.select_atoms("around "+str(r_c)+" (name ST* or name AT*)").ids)
    frame1 = frame1 + na
    #plot_group(frame1)
    #plt.title("Frame1")
     
    W.write(frame1) 
    #plot_group(frame1)
    #plot_group(frame1.select_atoms("name ST* or name AT*"))
    #plot_group(frame1.select_atoms("name O*"))
    #plt.show()
    #plot_group(frame1)
    #plot_group(frame1.atoms.select_atoms("name O*"))
    #plot_group(frame1.atoms.select_atoms("name ST* or name AT*"))
    #print(list(frame1.atoms.resnames))

    ###############FRAME 2, Half system adsorbed
    
    dim_x = u.trajectory[0].dimensions[0]
    dim_y = u.trajectory[1].dimensions[1]
    dim_z = u.trajectory[2].dimensions[2]
    cwls.atoms.positions = cw_orig_pos 
    
    sur_w = cwls.select_atoms("prop x > "+str(dim_x/2) )
    frame2 = cmls+sur_w
    ads_n = frame2.select_atoms("around "+str(r_c)+" (name ST* or name AT*)")
    print("2:")
    print(len(np.unique(ads_n.ids)))
    na = u.atoms - frame2
    na.atoms.positions = set_z_away(na)
    frame2 = frame2+na
    

    plot_group(frame2)
    #plt.show()
    W.write(frame2)
    #print(cwls.atoms.positions[2])
    #print(cwls.atoms.positions[2])
#############FRAME 3 all adsorb
    cwls.atoms.positions = cw_orig_pos 
    
    frame3 = cmls + cwls
    ads_n = frame3.select_atoms("(around "+str(r_c)+" (name ST* or name AT*)) and (not (name ST* or name AT* ))").n_atoms
    
    print("3:")
    print(ads_n)
    
    cwls.atoms.positions = cw_orig_pos 
    na = u.atoms - frame3
    na.atoms.positions = set_z_away(na)
    frame3 = frame3+na
    plot_group(frame3)
    #plt.show()
    W.write(frame3) 
    #FRAME 4  half adsorbed desorb
    cwls.atoms.positions = cw_orig_pos 


    frame4 = cmls + sur_w #  sur_w = cwls.select_atoms("prop x > "+str(dim_x/2) + " and prop y > "+str(dim_y/2))
    ads_n = frame4.select_atoms("(around "+str(r_c)+" (name ST* or name AT*)) and not (name ST* or name AT* )").n_atoms
    print("4:")
    print(ads_n)
    na = u.atoms - frame4
    na.atoms.positions = set_z_away(na)
    frame4 = frame4+na
    plot_group(frame4)
    #plt.show()
    W.write(frame4) 
    #Frame 5 Half adsorb, half adsorbed deorb
    cwls.atoms.positions =  cw_orig_pos

    frame5 = cmls + cwls.select_atoms("prop x < "+str(dim_x/2))
    ads_n = frame5.select_atoms("(around "+str(r_c)+" (name ST* or name AT*)) and not (name ST* or name AT* )").n_atoms
    print("5:")
    print(ads_n)
    
    na = u.atoms - frame5
    na.atoms.positions = set_z_away(na)
    frame5 = frame5+na
    W.write(frame5) 
    plot_group(frame5)
    #plt.show()
    #ALL DESORB
   
    cwls.atoms.positions = cw_orig_pos 
 
    frame6 = cmls
    ads_n = frame6.select_atoms("(around "+str(r_c)+" (name ST* or name AT*)) and not (name ST* or name AT* )").n_atoms
    print("6:")
    print(ads_n)
    na = u.atoms - frame6
    na.atoms.positions = set_z_away(na)
    frame6 = frame6+na
    plot_group(frame6)
    #plt.show()
    W.write(frame6)




