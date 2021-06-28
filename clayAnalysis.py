"""
TODO change print statments to logging
TODO allow specification of clay residues at instantiation?
"""
import MDAnalysis as mda
import nglview as nv
from nglview.datafiles import PDB, XTC
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def get_minima_coords(xs,ys):
        """
        Helper function that finds local minima of ys data and returns coordinates of these


        """
        
        wrap_len=len(ys)
        minima_coords = []
        for i in range(wrap_len):
            if (
             ((ys[i]==0 or (ys[i]< ys[(i+1)%wrap_len])) and
               ys[i]< ys[(i-1)%wrap_len]) 
             or 
             ((ys[i] < ys[(i+1)%wrap_len]) and
               ys[i]==0 or ys[i]< ys[(i-1)%wrap_len])): 
                minima_coords.append(xs[i])
        return minima_coords 



class ClayAnalysis:


    def __init__(self,u):
        self.universe = u
        self.box_dims = self.get_box_dim()#start stop etc should be defined



###Mostly stolen from dyndol (should import instead)########
    def get_box_dim(self,timestep=10, start=0,stop=-1):
    	'''
    	params: universe
    	params: timestep
    	return: box dimensions over time
    	See line 365 dynDol
    	'''
    	#logger.info("> getting simulation box dimensions...")
    	v = []
    	for ts in self.universe.trajectory[start:stop]:
    	    ptmp = self.universe.atoms.positions[:, 2]
    	    minpos = np.min(ptmp)
    	    maxpos = np.max(ptmp)
    	    dim = np.max(ptmp) - np.min(ptmp)
    	    v.append([ts.frame, minpos, maxpos, dim])

    	    #logger.debug(">> frame %s: z = %5.2f A..."%(ts.frame, dim))

    	box_dims = np.array(v)
    	box_dims[:, 0] *= timestep/1000.0

    	np.savetxt("bkp_box_dims.dat", box_dims) #Save dimensions to file
    	return box_dims


        	
    def get_partial_density(self, box_dims, sel, bins=100, start=0, stop=-1):
        '''        
        params: universe
        params: box_dims measuring box dimensions along z [frame, min, max, dim]
        params: selection name of residue(s) of interest
        params: number of bins
        params: start first frame to study
        params: stop last frame to study
        returns: density at each timepoint
        '''
        
        mol = self.universe.select_atoms(sel)

        res = [] #density collector
        cnt = 0
        for ts in self.universe.trajectory[start:stop]:

            curr_box = box_dims[cnt]
            binning = np.linspace(curr_box[1]-1, curr_box[2]+1, bins)
            
            result = np.histogram(mol.atoms.positions[:, 2], weights=mol.masses, bins=binning)    
            res.append(result)

            #logger.debug(">> frame %s density..."%ts.frame)
            cnt += 1

        return result





    def get_resnames_in_model(self):
        '''
        params: universe
        returns: unique residues in the universe
        '''
        
        all_labels = np.unique(self.universe.residues.resnames)
        #print(all_labels)
        all_select = []
        for l in all_labels:
            all_select.append("resname %s"%l)
        #print(all_labels)
        all_labels = np.concatenate((all_labels, ["system"]))
        all_select = np.concatenate((all_select, ["all"]))
        return all_select


###End of Mostly stolen from dyndol (should import instead)########


   





    def plot_average_density_t(self,res_density_w_time):#,z_max=133.69):
        """
        Plots the time averaged density distribution of all passed
        res_densitys from DynDen function output
        TODO generalise, Makes use of partial densities function of DynDen, make use of binning
        from r[0] info
        """
    
        t_steps = len(res_density_w_time[0])
        z_steps = len(res_density_w_time[0][0])
        #Bin overall size needs to be linked to the max size of system and number of points
    
        r_avgs = []
        for r in res_density_w_time:
    
            #Make use of existing bin data to plot
            r_avg_t = np.array(r[0])
    
            z_points = (r[1][1:] + r[1][:-1]) / 2
    
            for t in range(1,len(r)):
                r_avg_t = r_avg_t + r[1][t]
    
            plt.plot(z_points,r_avg_t/t_steps)




############Determining basal surfaces


    def generate_surface_groups(self,clay_resnames,box_dims=None):
        """
        params: clay_residues of interest to find surface of
        params: (optional) box dimension over time, default uses existing  
        returns: atomgroup of filtered passed clay residues that are at a basal surface

            For each residue passed, generate an atomgroup of those at the surface split into residues
            Uses the average density over time to calculate using minima of the edges of these
            
            TODO: would it be better these were just combined into one AtomGroup? Not a big change
            either way
            TODO: Also get interior surfaces
            TODO: Change over time might be an issue if clay not approx static
        """
        print("Determining basal surfaces of the clay:") 
        if(box_dims == None):
            box_dims = self.box_dims
        

        cly_comp = []
        density=[]
        lower_layers = []
        upper_layers = []
        
        for resname in clay_resnames:
            sel = "resname " +str(resname)
            cly_comp.append(self.universe.select_atoms(sel))
            density.append(self.get_partial_density(box_dims, sel))
    
            try:
                minima = get_minima_coords(density[-1][1],density[-1][0])
                
                #Want to find the z axis boundary of the lower clay basal surface
                lower_surface_bounds = minima[:2]
                #Want to find the z axis boundary of the other clay basal surface
                upper_surface_bounds = minima[-2:]
                
                
                #Use these found surface boundaries to create an atomgroup of upper and lower surface
                lower_layers.append(self.universe.select_atoms("(resname "+str(resname)
                                                 + ") and (prop z >= "+str(lower_surface_bounds[0])\
                                                 + ") and (prop z < "+str(lower_surface_bounds[1])\
                                                 +")"))
                
                upper_layers.append(self.universe.select_atoms("(resname "+str(resname)\
                                                 + ") and prop z >= "+str(upper_surface_bounds[0])\
                                                 +  " and prop z < "+str(upper_surface_bounds[1])))
                print(sel+" done")
            
            except IndexError: #If density is entirely 0
                print("Error looking at surfaces of  resname " + resname + " (might be caused by it having no density or obvious minimums?)")
    
        
        return [lower_layers,upper_layers]
        
        

###############################Adsorption analysis


    def find_adsorbed(self,surface_ids,adsorbants,r_c):
        """
        params: Indices of atoms at the surface of clay (AtomGroup.indices)
        params: Name of possible adsorbants to the surface of clay that are of interest
        params: Cutoff distance at which a particle is considered to be adsorbed
        returns: adsorbed atoms of interest at current time step
    
        Create a list of search strings that filters atoms in a certain distance to the
        surface and of a specific type of interest
        """
    
        search_strs = ""
    
        for i in surface_ids:
            for j in adsorbants:
                search_strs+= "(resname " + str(j) + " and around "+ str(r_c) + " index "+str(i)+") or "
                #Might be quicker to do as two lists
        #print(search_strs[:-3])
        return self.universe.select_atoms(search_strs[:-3])
    
    def find_adsorbed_as_agroup(self,surface_ids,adsorbants,r_c):
        """
        params: index of atoms at the surface of clay
        params: resnames of possible adsorbants to the surface of clay that are of interest
        params: Cutoff distance at which a particle is considered to be adsorbed
        returns: list of atomGroups of adsorbed ions (of type given in adsorbant) to given surface_ids atoms
    
        Creates a 2D list of atomgroups, which links a surface atom to a number of adsorbant ions
        """
    
        adsorbed_to_surface = []
    
        for i in range(len(surface_ids)):
            search_strs= ""
            for p_ads in adsorbants:
                search_strs+= "(resname " + str(p_ads) + " and around "+ str(r_c) + " index "+str(surface_ids[i])+") or "
                #Might be quicker to do as two lists
            adsorbed_to_surface.append([self.universe.select_atoms("index "+str(surface_ids[i])),u.select_atoms(search_strs[:-3])])
    
        return adsorbed_to_surface
