"""
TODO change print statments to logging
TODO allow specification of clay residues at instantiation?
TODO Get rid of box_dims for timestep.dimensions
"""
import MDAnalysis as mda
import nglview as nv
from nglview.datafiles import PDB, XTC
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from MDAnalysis import transformations
from MDAnalysis import analysis
import scipy.stats as stats
from scipy.signal import find_peaks



def get_minima_coords(xs,ys):
        """
        Helper function that finds local minima of ys data and returns coordinates of these


        """
        #ys = ys*-1
        #x_peak_ind = find_peaks(ys)
        #print("peak finding")
        #print(x_peak_ind)
        #print(xs[x_peak_ind])
        #return xs[x_peak_ind]
         

        
        

        wrap_len=len(ys)
        minima_coords = []
        for i in range(wrap_len):
            if(ys[i]<ys[(i+1)%wrap_len] and ys[i] < ys[(i-1)%wrap_len]):
                    minima_coords.append(xs[i])

            elif(ys[i]==0 and (ys[(i+1)%wrap_len]>0 or  ys[(i-1)%wrap_len]>0)):
                    minima_coords.append(xs[i])
        return minima_coords 


def shift_clay_to_center(ts):
    """
    Translate all coordinates to have clay centered.

    Might be better to pass clay atom selection? simpler version for now
    """
    clay = u.select_atoms("resname NON*")
    
    
    ts.positions = ts.positions - clay.center_of_mass()
    return ts


class ClayAnalysis:


    def __init__(self,u):
        self.universe = u
        self.box_dims = self.get_box_dim()#start stop etc should be defined

    def shift_clay_to_origin(self, u):
        u.trajectory.add_transformations(*[shift_clay_to_center,mda.transformations.wrap(u.atoms)])

###Mostly stolen from dynden (should import instead)########
    def get_box_dim(self,timestep=10, start=0,stop=-1):
    	'''
    	params: universe
    	params: timestep
    	return: box dimensions over time
    	See line 365 dynden
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
        #bins = int(np.sqrt(mol.n_atoms))
        #print(int(bins))
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


###End of Mostly stolen from dynden (should import instead)########


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
        ids = ["Non OW", "STs", "ATs"]
        r_avgs = []
        i = 0 
        for r in res_density_w_time:
            
            #print(r) 
            #Make use of existing bin data to plot
            r_avg_t = np.array(r[0])
             
    
            for t in range(1,len(r)):
                r_avg_t = r_avg_t + r[0][t]
   


            zpoints = (r[1][1:] + r[1][:-1]) / 2 #Use midpoint of bins
            #print(zpoints)

            #print((len(zpoints),len(r[1])))
            plt.plot(zpoints,r[0])
                        



############Determining basal surfaces

    def generate_surface_group(self,surf_type="waters",box_dims=None):
        """ 
        #params: clay components of interest e.g. surface waters or surfce minerals 
        #params: (optional) box dimension over time, default uses existing  
        #returns: atomgroup of filtered clay components that are at a basal surface

            
        #    TODO: Change over time might be an issue if clay not approx static
        """
        print("Determining basal surfaces of the clay:") 
        if(box_dims == None):
            box_dims = self.box_dims
        

        cly_comp = []
        density=[]
        lower_layers = []
        upper_layers = []
        
        if(surf_type == "mineral"): 
                atom_sel =["name ST* or name AT*"]

                den = self.get_partial_density(box_dims, "resname NON* and (name ST* or name AT*)")#or name AT*")
                peakNum = [5,7]
        else:
                atom_sel = ["resname NON* and name O*"]#["name O[!W]* and name O[!T]*"]
                den = self.get_partial_density(box_dims," resname NON* and name O*") 
                peakNum = [6,8]
        #hist_boxes_to_ midpoints
        z_points = (den[1][1:] + den[1][:-1]) / 2
        minima = get_minima_coords(z_points,den[0])
        #Want to find the z axis boundary of the lower clay basal surface

        lower_surface_bounds = minima[peakNum[0]-1:peakNum[0]+1]
        

        print(lower_surface_bounds)
        #Want to find the z axis boundary of the other clay basal surface
        upper_surface_bounds = minima[peakNum[1]-1:peakNum[1]+1]
        print(upper_surface_bounds)
        
        for sel in atom_sel:
            cly_comp.append(self.universe.select_atoms(sel))
            
            try:
                
                
                #Use these found surface boundaries to create an atomgroup of upper and lower surface
                lower_layers.append(self.universe.select_atoms("("+sel
                                                 + ") and (prop z >= "+str(lower_surface_bounds[0])\
                                                 + ") and (prop z < "+str(lower_surface_bounds[1])\
                                                 +")"))
                
                upper_layers.append(self.universe.select_atoms("("+sel
                                                 + ") and prop z >= "+str(upper_surface_bounds[0])\
                                                 +  " and prop z < "+str(upper_surface_bounds[1])))
            
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
