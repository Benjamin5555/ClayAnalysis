"""
TODO:
    *change print statments to logging
    *allow specification of clay residues at instantiation?
    *Get rid of box_dims for timestep.dimensions
    *make transformation of clay to origin automatic
    *get_partial_density function currently only returns a single partial density (last frame) need
     to change it's return to res and then update rest of code to deal with this change 

"""
import MDAnalysis as mda
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from MDAnalysis import transformations
from MDAnalysis import analysis
import scipy.stats as stats
#from scipy.signal import find_peaks
import sys
import csv

def hist(x_data,bins=None,x_title=None,y_title=None,y_labels=None,title=None):
    """
    Plot a histogram of some passed data

    Wrapper fo matplotlib histogram plotting

    Args:
        x_data: Data to be histogramed
        bins:   (Optional) number of bins to use
        x_title:(Optional) Label for x_axis of plot
        y_title:(Optional) Label for y_axis of plot
        title:  (Optional) Plot title
    Returns:
        Histogram values and bins
    Raises:
        ZeroDivisionError, AssertionError, & ValueError.

    Examples:

    """

    if(bins==None):
        bins = len(x_data)/10

    ns,bins,*_ = plt.hist(x_data,bins, alpha=0.75)    

    # Set axes, titles and legend
    plt.title(title)
    plt.ylabel(y_title)
    plt.xlabel(x_title)
    plt.legend()
    #plt.show()
    if(title != None):
        plt.savefig(str(title)+".png")
    return ns,bins

def plot_group(grp,fig=None,ax=None,label=None):
    """
    Plot positions of atoms in a passed group

    Creates simple 3D plot of the positions of atoms in the passed group and returns a figure 
    and axes object on which more groups can be plotted

    Args:
        grp:    Group to plot positions of
        fig:    (Optional) Figure to add plots to
        ax:     (Optional) Axes to add plots to
        label:  (Optional) Label to add to plot for passed group
    Returns:
        Figure and axes containing the plot of atomgroup positions

    Raises:
        ZeroDivisionError, AssertionError, & ValueError.

    Examples:
        >>> atom_group = u.atoms.select_atoms("resname NON*")
        >>> plot_group(atom_group)
        >>> plt.show()
    """
    if fig ==None:
        fig = plt.figure()
        ax = Axes3D(fig, zlabel="z")

    upper_x = np.transpose(grp.positions)[0]
    upper_y = np.transpose(grp.positions)[1]
    upper_z = np.transpose(grp.positions)[2]

    ax.scatter(upper_x,upper_y,upper_z,label=label)
    return fig, ax

def get_minima_coords(xs,ys):
    """
    Finds local minima of ys data and returns xs coordinates of these

    Simple minima finder based on minimums or starts of zeroed data

    Args:
        xs:   X coordinates
        ys:   related Y coords to find minimums of


    Returns:
        List of positions in xs that relate to minimums in ys 

    Raises:
        ZeroDivisionError, AssertionError, & ValueError.

    Examples:
        >>> xs = [0,1,2,3,4,5,6,7]
        >>> ys = [0,0,0,5,4,1,4,5]
        >>> get_minima_coords(xs,ys)
        [2,5]
    
    """
    wrap_len=len(ys)
    minima_coords = []
    for i in range(wrap_len):
        if(ys[i]<ys[(i+1)%wrap_len] and ys[i] < ys[(i-1)%wrap_len]):
                minima_coords.append(xs[i])

        elif(ys[i]==0 and (ys[(i+1)%wrap_len]>0 or  ys[(i-1)%wrap_len]>0)):
                minima_coords.append(xs[i])

    return minima_coords 




class ClayAnalysis:


    def __init__(self,u):
        """
        Initalise the clay analysis object, by linking a universe object to it 

        Required to link a universe that is under analysis as well as provide preprocessing to 
        the system, e.g. by moving the clay to the origin of the simulation box.

        Class provides some functionality to provide information about adsorption to the surface
        of a clay e.g. determining a list of adsorption times and statistics
        
        Args:
            u:      Universe object

        Raises:
            ZeroDivisionError, AssertionError, & ValueError.

        Examples:
            >>> u = mda.Universe("topo.tpr","traj.trr")
            >>> cal = ClayAnalysis.clayAnalysis(u)
        """
        #TODO Might be better if this was setup to just extend the universe object 

        self.universe = u
        self.box_dims = self.get_box_dim()#start stop etc should be defined
        self.clay = self.universe.select_atoms("resname NON*")
        self.move_clay_to_origin()
        
    def combine_atomgroups(self,list_of_groups):
        """
        Converts a list of atomgroups into a single atomgroup
        
        Args:
            list_of_groups:      Python list made up of atomgroup objects

        Returns:
            Single atomgroup containing all atoms in passed list

        Examples:
            >>> a = u.select_atoms("name AtomA")
            >>> b = u.select_atoms("name AtomB")
            >>> cal = ClayAnalysis.clayAnalysis(u)
            >>> combined_group = cal.combine_atomgroups([a,b])
            >>> combined_group.names
            ["AtomA","AtomB"]
        """

        
        ag = mda.AtomGroup([],self.universe)#list_of_groups[0]
        for g in list_of_groups:
            ag = ag + g
        return ag


    def __shift_clay_to_origin(self,ts):
        """
        Transformation function that moves all coordinates to have clay centered at the origin 

        Should be used only with .add_transformation function
        """
        def wrapped(ts):
            ts.positions = ts.positions - self.clay.center_of_mass()
            return ts
        return wrapped(ts)


###Mostly stolen from dynden (should import instead)########
    def get_box_dim(self,timestep=10, start=0,stop=-1):
        """
        Creates a list of simulation box dimensions with time

        Required for use with the partial density function that was taken from DynDen 
        (https://github.com/punkpony/DynDen.git)

        Args:
            timestep:  (Optional) Timestep size, default 10     
            start:     (Optional) starting timstep to analyse, default 0
            stop:      (Optional) Final timestep to analyse, default -1 (final step)


        Returns:
            Simulation box dimensions over time
        
        Raises:
            ZeroDivisionError, AssertionError, & ValueError.

        Examples:
            >>> cal = ClayAnalysis.clayAnalysis(u)
            >>> box_dims_w_time = cal.get_box_dim()

        """

        

        #TODO: This can probably be simplified into a loop over trajectroy creating a list of u.dimensions objects
        v = []
        if (self.universe.trajectory.n_frames>1):
            for ts in self.universe.trajectory[start:stop]:
                ptmp = self.universe.atoms.positions[:, 2]
                minpos = np.min(ptmp)
                maxpos = np.max(ptmp)
                dim = np.max(ptmp) - np.min(ptmp)
                v.append([ts.frame, minpos, maxpos, dim])


        else:
            ptmp = self.universe.atoms.positions[:, 2]
            minpos = np.min(ptmp)
            maxpos = np.max(ptmp)
            dim = np.max(ptmp) - np.min(ptmp)
            v.append([self.universe.trajectory.frame, minpos, maxpos, dim])

        box_dims = np.array(v)
        box_dims[:, 0] *= timestep/1000.0

                #np.savetxt("bkp_box_dims.dat", box_dims) #Save dimensions to file
        return box_dims


                
    def get_partial_density(self, box_dims, sel, bins=100, start=0, stop=-1,density=False):
        """
            Creates a list of simulation box dimensions with time

            Required for use with the partial density function that was taken from DynDen 
            (https://github.com/punkpony/DynDen.git)

            Args:
                box_dims: Simulation box dimensions for each time step 
                sel:      Selection string for atoms using MDAnalysis selection string 
                bins:     (Optional) Number of bins to use, default 100
                start:    (Optional) starting timstep to analyse, default 0
                stop:     (Optional) Final timestep to analyse, default -1 (final step)


            Returns:
                Density profile along z-axis (assumed normal to clay) for each timepoint
                in the form of a list of tuples of (histogramValues,binsizes) for each timestep 
                
            Raises:
                ZeroDivisionError, AssertionError, & ValueError.

            Examples:

        """
        print("Calculating partial densitiy for selection " +str(sel))
        mol = self.universe.select_atoms(sel)
        res = [] #density collector
        #print(mol)   
        cnt = 0
        if(self.universe.trajectory.n_frames>1):
                for ts in self.universe.trajectory[start:stop]:

                    curr_box = box_dims[cnt]
                    binning = np.linspace(curr_box[1]-1, curr_box[2]+1, bins)
                   
                    result = np.histogram(mol.atoms.positions[:, 2], weights=mol.masses, bins=binning)

                    res.append(result)

                    #logger.debug(">> frame %s density..."%ts.frame)
                    cnt += 1
        else:
                curr_box = box_dims[0]
                binning = np.linspace(curr_box[1]-1, curr_box[2]+1, bins)
                    
                result = np.histogram(mol.atoms.positions[:, 2], weights=mol.masses, bins=binning)
                res.append(result)
        return res



    #def get_resnames_in_model(self):
    #    """
    #        Simply returns a selection string for all resn

    #        Required for use with the partial density function that was taken from DynDen 
    #        (https://github.com/punkpony/DynDen.git)

    #        Args:
    #            box_dims: Simulation box dimensions for each time step 
    #            sel:      Selection string for atoms using MDAnalysis selection string 
    #            bins:     (Optional) Number of bins to use, default 100
    #            start:     (Optional) starting timstep to analyse, default 0
    #            stop:      (Optional) Final timestep to analyse, default -1 (final step)


    #        Returns:
    #            Density profile along z-axis (assumed normal to clay) for each timepoint
    #            in the form of a list of tuples of (histogramValues,binsizes) for each timestep 
    #            
    #        Raises:
    #            ZeroDivisionError, AssertionError, & ValueError.

    #        Examples:

    #    """
    #    
    #    all_labels = np.unique(self.universe.residues.resnames)
    #    for l in all_labels:
    #        all_select.append("resname %s"%l)
    #    all_labels = np.concatenate((all_labels, ["system"]))
    #    all_select = np.concatenate((all_select, ["all"]))
    #    return all_select


###End of Mostly stolen from dynden (should import instead)########
                        
    def plot_density_w_time(self,selection,label=""):
        """
            Creates a matplotlib object of the density profile averaged over time for a selection

            Given an atomgroup object, plots the time averaged mass density of these atoms along the
            z-axis 

            Args:
                sel:      Atomgroup selection of atoms to plot the time averaged mass density of 
                label:    Label for the key of the plot object for this selection

            Examples:

        """

        den = self.get_partial_density(self.box_dims," or ".join(["name "+s for s in np.unique(selection.names)]))
        z_points = (den[0][1][1:] + den[0][1][:-1]) / 2
        for t_step in range(len(den)):
            sum_den = den[t_step][0]
        avg_den = sum_den/len(den[0])
        plt.plot(z_points,avg_den,label=label)
        plt.savefig(str(label)+" density.png")

############Determining basal surfaces

    def move_clay_to_origin(self):
        """
            Shift entire clay such that its center of mass is at 0 along the z-axis

            Is called upon setup of a clay analysis object and so is not really needed to be called 
            directly


            Examples:

        """
        print("Shifting clay center of mass to be at origin")
        try:
            self.universe.trajectory.add_transformations(*[self.__shift_clay_to_origin,mda.transformations.wrap(self.universe.atoms)])
        except:
            print("Unable to shift the clay to the origin, perhaps already shifted?") 


    def generate_surface_group(self,surf_type="waters",box_dims=None):
        """
            Creates a list of two atomgroups, which represent the 'top' and 'bottom' basal surfaces 

            Returns two atomgroups, the lower and upper surface respecively.
            NOTE:Makes use of assumption of 3 layered clay and uses the mass density along z to work
            out which atoms are at surface. As such if the density profile is significantly 
            different in terms of number of peaks, this code will identify the incorrect groups as 
            surface atoms

            Args:
                surf_type: (optional) clay components of interest options include waters or mineral,
                            default waters (split due to differing mass density profile) 

                box_dims: (optional) box dimension over time, default uses existing  
        

            Returns:
                List of 2 atomgroups of clay components that are at a basal surface, first element 
                being the lower surface, the second the upper surface

                
            Raises:
                ZeroDivisionError, AssertionError, & ValueError.

            Examples:

        """ 

            
        #TODO: Change over time might be an issue if clay not approx static
        print("Determining basal surfaces of the clay:") 

        try:
            self.move_clay_to_origin()
        except ValueError:
            print("Clay already centered to origin")
        if(box_dims == None):
            box_dims = self.box_dims
        
        cly_comp = []
        density=[]
        lower_layers = []
        upper_layers = []
        
        if(surf_type == "mineral"): 
                atom_sel =["name ST* or name AT*"] #Should just be not certain types
                den = self.get_partial_density(box_dims, "resname NON* and (name ST* or name AT*)")#or name AT*")
                peakNum = [5,7]
        else:
                atom_sel = ["resname NON* and name O*"]#["name O[!W]* and name O[!T]*"]
                den = self.get_partial_density(box_dims," resname NON* and name O*") 
                #self.plot_density_w_time(self.universe.atoms.select_atoms("resname NON* and name O*"))


                peakNum = [5,7]

        #Makes use of initial state only
        z_points = (den[0][1][1:] + den[0][1][:-1]) / 2
        minima = get_minima_coords(z_points,den[0][0])
        #Want to find the z axis boundary of the lower clay basal surface

        lower_surface_bounds = minima[peakNum[0]-1:peakNum[0]+1]
        
        #Want to find the z axis boundary of the other clay basal surface
        upper_surface_bounds = minima[peakNum[1]-1:peakNum[1]+1]
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
                print("Error looking at surfaces of selection " + sel + " (might be caused by it having no density or obvious minimums?)")
               
        
        return [lower_layers,upper_layers]





###############################Adsorption analysis


    def find_adsorbed(self,surface_ids,adsorbants_resnames,r_c_upper,r_c_lower=0):#,r_c=10):
        """
            Return a dictionary of surface ions and the things adsorbed to them at the current step

            By looping over every surface atoms, generates a dictionary of surface atoms (atomgroup 
            with single surface atoms within) with something adsorbed to them to each cation that is
            adsorbed using the radial disance between the surface atom and cation as a cutoff of 
            adsorbed versus non adsorbed


          
            Args:
                surface_ids: index of atoms at the surface of clay
                adsorbants_resnames: resnames of possible adsorbants to the surface of clay that are
                 of interest
                r_c_upper: Cutoff distance at which a particle is considered to be adsorbed
                r_c_lower: Distance along z above which a particle is considered to be adsorbed

            Returns:
                Dictionary of atomGroups of adsorbed ions (of type given in adsorbant) to given 
                surface atoms
                i.e. {Surface_atom (Atomgroup containing 1 atom): [Adsorbed_cations (Atomgroup)]} 

                
            Raises:
                ZeroDivisionError, AssertionError, & ValueError.

            Examples:

        """ 
        adsorbed_to_surf = {}
        for i in range(len(surface_ids)):
            
            search_strs= ""
            adsorbed = self.universe.select_atoms("")
            
            for p_ad in adsorbants_resnames: #Don't think need to do the around thing every time here?
                search_strs+= "(resname " + str(p_ad) + " and sphlayer "+ str(r_c_lower) + " "+str(r_c_upper) +" index "+str(surface_ids[i])+") or "
                #search_strs+= "(resname " + str(p_ad) + " and around "+ str(r_c) + " index "+str(surface_ids[i])+") or "
                #Might be quicker to do as two lists

            surf_id = self.universe.select_atoms("index "+str(surface_ids[i]))
            adsorbed = self.universe.select_atoms(search_strs[:-3])
            #print(adsorbed)
            if (not adsorbed.n_atoms == 0):  
                
                #If not adsorbed to anything, don't bother adding to dict
                adsorbed_to_surf[surf_id]= adsorbed
        return adsorbed_to_surf

    def find_adsorption_times(self,lower_surf_grp, upper_surf_grp, adsorbant_resname, r_c_upper=0, r_c_lower=0, start=0, stop=-1):
        """
            Determines the times overwhich an ions are adsorbed to a surface by looking at z position

            Simply uses z position of ions relative to the surface of the clay as a definition for 
            adsorption. i.e. if an ion has position along the z-axis a distance between r_c_lower 
            and r_c_upper to the average position of surface clay atoms

            ##REPLACE BELOW WITH AN IMAGE 
            upper surface
            -------------
            r_c_lower
            -------------
            adsorbed ions
            -------------
            r_c_upper
            -------------
            Bulk Solution
            -------------
            r_c_upper
            -------------
            adsorbed ions
            -------------
            r_c_lower
            -------------
            lower surface
            
            NOTE: If an ion is 'already adsorbed' at the start of the simulation it will be 
                  considered as though it adsorbs at starting time step. As such ions should start 
                  off far from the clay surface (achiveable via the reinsert_ions.py script in 
                  https://github.com/Benjamin5555/Clay-Project)
            NOTE: If an ion is still adsorbed at the end of the simulation its time will not be 
                  output 

            Args:
                lower_surf_grp: 'bottom' clay surface atomgroup
                upper_surf_grp: 'top' clay surface atomgroup
                adsorbants_resnames: resnames of possible adsorbants to the surface of clay that are
                                     of interest
                r_c_upper: Cutoff distance at which a particle is considered to be adsorbed
                r_c_lower: Distance along z above which a particle is considered to be adsorbed
                start:     (Optional) starting timstep to analyse, default 0
                stop:      (Optional) Final timestep to analyse, default -1 (final step)

            Returns:
                Statistics on adsorption and a list of adsorption times collected over the 
                trajectory for atoms that are within the distance boundary  
                
                Statistics have form Step,Number of  adsorption events, num desorption events, 
                number continuing to be adsorbed, number adsorbed total at ts, total_adsorption 
                events over whole trajectory
 
            Raises:
                ZeroDivisionError, AssertionError, & ValueError.

            Examples:
                >>> cal = clayAnalysis.ClayAnalysis(u)
                >>> surfaces = cal.generate_surface_group("mineral")
                >>> lower_surf_grp = cal.combine_atomgroups(surfaces[0])
                >>> upper_surf_grp = cal.combine_atomgroups(surfaces[1])
                >>> times,stats = cal.find_adsorbant_times(lower_surf_grp,upper_surf_grp,"Cs",4.8,0)
        """ 
        adsorbant_resname = np.unique(adsorbant_resname)
        self.universe.trajectory[0]
        current_adsorbed = self.universe.select_atoms("")
        historic_adsorbed = {}
        times = [] 
        stats = []
        tot_ads = 0

        avg_surf_b = np.average(upper_surf_grp.positions.T[2])
        avg_surf_a = np.average(lower_surf_grp.positions.T[2])
        surf_b_ub = avg_surf_b-r_c_lower
        surf_b_lb = avg_surf_b-r_c_upper 
        surf_a_ub = avg_surf_a+r_c_upper 
        surf_a_lb = avg_surf_a+r_c_lower

        with open(str(adsorbant_resname[0])+'stats.csv', 'w', newline='') as csvfile:
                statWriter = csv.writer(csvfile, delimiter=',')

                for ts in self.universe.trajectory:
                    print(adsorbant_resname[0])

                    #counters 
                    c_ad = 0  
                    c_ct = 0 
                    c_ds = 0
                    
                    current_adsorbed_ls = self.universe.select_atoms("resname " + str(adsorbant_resname[0])+" and prop z > "+str(surf_a_lb)+" and prop z < "+str(surf_a_ub))
                    current_adsorbed_us = self.universe.select_atoms("resname " + str(adsorbant_resname[0])+" and prop z > "+str(surf_b_lb)+" and prop z < "+str(surf_b_ub))
                    
                    
                    current_adsorbed = current_adsorbed_ls + current_adsorbed_us 
                    hist_ads_grp = self.combine_atomgroups(historic_adsorbed.keys())
                    desorbed = hist_ads_grp - current_adsorbed
                    adsorbed = current_adsorbed - hist_ads_grp
                    cont_ads = current_adsorbed & hist_ads_grp
                    
                    for ads_id in adsorbed:
                        historic_adsorbed[ads_id] = 1 # Add a time adsorbed for
                        c_ad = c_ad + 1
    
                    for ct_id in cont_ads:
                        #NOTE: If instead you recorded the timestep at which this ion was adsorbed and used the difference to the step it desorbs, you would not require this incrementing step over all ions
                        historic_adsorbed[ct_id] += 1 
                        c_ct = c_ct + 1

                    for desorbed_id in desorbed:
                        times.append(historic_adsorbed.pop(desorbed_id))
                        c_ds = c_ds + 1
                    
                    tot_ads = tot_ads + c_ad
                    stats.append([c_ad,c_ds,c_ct,c_ad+c_ct,tot_ads]) 
                    print("STATS OUTPUT"+str(self.universe.trajectory.ts))
                    print([self.universe.trajectory.time]+stats[-1])
                    statWriter.writerow(stats[-1])

        np.savetxt(str(adsorbant_resname[0])+"OuputTimes.txt",times)
        return times, stats


    def find_adsorption_times_c(self,surface_ids,adsorbant_resname,r_c_upper=0,r_c_lower=0,start=0,stop=-1):
        """
            Determines the times overwhich ions are adsorbed to a surface by looking at radial 
            distance of all surface atoms

            Uses a similar process to find_adsorpion_times in terms of looking at changes in groups 
            but is able to keep track of specific surface atom to specific ion relationships and so
            could be built upon to gain furher info on specific adsorption events, should it ever be 
            needed


                       
            NOTE: If an ion is 'already adsorbed' at the start of the simulation it will be 
                  considered as though it adsorbs at starting time step. As such ions should start 
                  off far from the clay surface (achiveable via the reinsert_ions.py script in 
                  https://github.com/Benjamin5555/Clay-Project)
            NOTE: If an ion is still adsorbed at the end of the simulation its time will not be 
                  output 

            Args:
                surface_ids: atomgroup.indices property of the surface group of interest to look at 
                             adsorbance to 
                lower_surf_grp: 'bottom' clay surface atomgroup
                upper_surf_grp: 'top' clay surface atomgroup
                adsorbants_resnames: resnames of possible adsorbants to the surface of clay that are
                                     of interest
                r_c_upper: Cutoff distance at which a particle is considered to be adsorbed
                r_c_lower: Distance along z above which a particle is considered to be adsorbed
                start:     (Optional) starting timstep to analyse, default 0
                stop:      (Optional) Final timestep to analyse, default -1 (final step)

            Returns:
                Statistics on adsorption and a list of adsorption times collected over the 
                trajectory for atoms that are within the distance boundary  
                
                Statistics have form Step,Number of  adsorption events, num desorption events, 
                number continuing to be adsorbed, number adsorbed total at ts, total_adsorption 
                events over whole trajectory
 
            Raises:
                ZeroDivisionError, AssertionError, & ValueError.

            Examples:
                >>> cal = clayAnalysis.ClayAnalysis(u)
                >>> surfaces = cal.generate_surface_group("mineral")
                >>> lower_surf_grp = cal.combine_atomgroups(surfaces[0])
                >>> upper_surf_grp = cal.combine_atomgroups(surfaces[1])
                >>> times,stats = cal.find_adsorbant_times(lower_surf_grp,upper_surf_grp,"Cs",4.8,0)
                        

        """
        #TODO: Allow different ions to have different adsorption cutoff distances
        #TODO: Case where ions adsorbed in first time step 
        #        ->considered to be adsorbed as though just adsorbed to surface
        #TODO: Case where ions still adsorbed to the surface at the end -> times not recorded
        times = [] 
        
        prev_ads_record ={}  
        # {surface_atom: {adsorbant:time adsorbed}}
        # Record of surface atoms and their adsorbants as well as time adsorbed
        
        stats = []
        counters = [0,0,0,0,0]
        tot_ads = 0  
        
        print("Step,Num ads, num ds, num continue, num adsorbed tot at ts,total_adsorption events")
        with open(adsorbant_resname+'stats.csv', 'w', newline='') as csvfile:
                statWriter = csv.writer(csvfile, delimiter=',')
 
                for ts in self.universe.trajectory:
                        #counters 
                        c_ad = 0  
                        c_ct = 0 
                        c_ds = 0
                        
                        #Get surface atoms and attached ions  
                        # {surface_atom: adsorbant}
                        currently_ads_dict = self.find_adsorbed(surface_ids,adsorbant_resname,r_c_upper,r_c_lower)
                        print(currently_ads_dict) 
                        #If a surface ion is no longer sorbed to anything, we remove it from the list of stored ions
                        #TODO Combine into below loop?
                        rm =   prev_ads_record.keys() - currently_ads_dict.keys()
                        for rem_surf_at in rm:
                            removed = prev_ads_record.pop(rem_surf_at)
                            for rem_ion in removed:
                                times.append(removed[rem_ion])
                                
                                c_ds = c_ds +1 
                       
                        #Look at each surface ion which has something adsorbed to it and compare to historic 
                        #record of adsorpion info (no interdependency other than stats)
                        for surf_atm in currently_ads_dict.keys():
                            
                            #if surface atom not previously adsorbed to then will need to be added to record
                            #dictionary
                            if(not surf_atm in prev_ads_record.keys()):
                                prev_ads_record[surf_atm] = {}
                              
                            #Combine together list of atoms adsorbed to the surface atom into a single 
                            #atomgroup (to perform group operations)
                            prev_adsorbed = self.combine_atomgroups(list(prev_ads_record[surf_atm].keys()))
                            prev_ads_record_at_surf_atm = prev_ads_record[surf_atm]

                        
                            currently_adsorbed = currently_ads_dict[surf_atm] 
                                            
                            #Looking at in current but not historic (newly adsorbed)
                            newly_adsorbed = currently_adsorbed - prev_adsorbed 
                            c_ad= c_ad + self._update_record_w_newly_adsorbed(newly_adsorbed,prev_ads_record_at_surf_atm )
                                
                            #In both historic and current (continue to be adsorbed)
                            continue_to_adsorb = currently_adsorbed & prev_adsorbed 
                            #NOTE: If instead you recorded the timestep at which this ion was adsorbed and used the difference to the step it desorbs, you would not require this incrementing step over all ions

                            c_ct = c_ct + self._update_record_w_continuing_adsorbed(continue_to_adsorb,prev_ads_record_at_surf_atm )
                            
                            #Only in historic (desorbed this step) 
                            newly_desorbed = prev_adsorbed - currently_adsorbed 
                            c_ds = c_ds + self._update_record_w_newly_desorbed(newly_desorbed,prev_ads_record_at_surf_atm,times) 
                            
                        tot_ads = tot_ads + c_ad
                        stats.append([c_ad,c_ds,c_ct,c_ad+c_ct,tot_ads]) 
                        
                        print([self.universe.trajectory.ts,self.universe.trajectory.time]+stats[-1])
                                        
                        statWriter.writerow(stats[-1])
        np.savetxt("ouputTimes.txt",times)
        return times, stats

    def _update_record_w_newly_adsorbed(self,newly_adsorbed,prev_ads_record_at_surf_atm):
        """
            Update 'historic record' of what is adsorbed at the previous timestep to add in any new 
            ions that adsorb at that step
        """

        c_ad = 0
        for ads in newly_adsorbed:
            c_ad = c_ad +1
            prev_ads_record_at_surf_atm[ads]=1 #Assumes this to be starting point +-dt error
        return c_ad
    
    def _update_record_w_continuing_adsorbed(self,continue_to_adsorb,prev_ads_record_at_surf_atm):
        """
            Update the record of ions that are continuing to be adsorbed to a surface ion by 
            increment the time for which they have been adsorbed
        """
        c_ct =0 
        for inc in continue_to_adsorb: 
            time = prev_ads_record_at_surf_atm.pop(inc)
            prev_ads_record_at_surf_atm[inc] = time+1#self.universe.dt
            c_ct = c_ct + 1
        return c_ct            

    def _update_record_w_newly_desorbed(self,newly_desorbed,prev_ads_record_at_surf_atm,times):
        """
            Update record for surface ions that were adsorbed to multiple ions but in the current 
            time step have had one or more (but not all) of these desorb. Recording time overwhich 
            they were adsorbed in # time steps
        """
        c_ds = 0 
        for dsb in newly_desorbed:
            times.append(prev_ads_record_at_surf_atm.pop(dsb)) #Assumes this to be ending point +-dt error
            c_ds = c_ds+1
        return c_ds 

    def get_num_adsorbed_at_current_time(self, adsorbed_dict):
        """
            Return total number of ions adsorbed to surface in a time step
        """
        count  = 0
        for surf_atm in adsorbed_dict.keys():
            for ads in adsorbed_dict[surf_atm]:
                count = count+1      
        return count  

#####################################CHARGE THROUGH A SURFACE

    def plot_surface_charge(self,selection,bins=100,interpolate=True):
        """
            params: Atomgroup of interest 
            params: (optional) bin sizes
            params: (optional) bool smoothed 
            
            Simple charge surface plotting via a 2d histogram
            NOTE: Not useful should just use VMD instead
        """

        Xs = selection.atoms.positions[:,0]
        Ys = selection.atoms.positions[:,1]
        Zs = selection.atoms.positions[:,2]
        fig, ax = plt.subplots()

        im = plt.hist2d(Xs,Ys,weights=selection.atoms.charges,bins=200)
        if(interpolate):
            ax.imshow(im[0], interpolation = "gaussian",extent=[0,self.universe.dimensions[0],0,self.universe.dimensions[1]])

        fig.colorbar(im[3])

        plt.show()




if __name__ == "__main__":
    if len(sys.argv) == 0:
        print("usage:")
        print("topology file, trajectory file, mode, mode specific options")
        print("Modes:")
        print("0 : Plot densities, mode specific options: selection strings e.g. 'resname NON*,resname Cs'")
    topfile = sys.argv[1] 
    trajfile = sys.argv[2]
    mode = sys.argv[3]
    if(topfile[-3:] =="top"):
        u = mda.Universe(topfile,trajfile,topology_format='ITP')
    else:
        u = mda.Universe(topfile,trajfile)
    #u = mda.Universe("TestFiles/test_sys.pdb","TestFiles/test_sys.pdb")
    cal = ClayAnalysis(u)
    
    if mode == "0" :
        #Density analysis
        sel_strs = sys.argv[4:]
        for sel_str in sel_strs:
            selection = u.atoms.select_atoms(sel_str)
            print(np.unique(selection.resnames))
            print(selection.n_atoms)
            cal.plot_density_w_time(selection,label=sel_str)
        plt.legend()
        
        plt.show()

    if mode == "2":
        print("Adsorption times analysis")
        #cal.move_clay_to_origin()
        r_c_lower = float(sys.argv[4])
        r_c_upper = float(sys.argv[5])
        ads_sel_str = sys.argv[6]
        ads = u.atoms.select_atoms(ads_sel_str)
        lower_layers,upper_layers = cal.generate_surface_group("mineral")
        
        lower_surf_grp = cal.combine_atomgroups(lower_layers)
        upper_surf_grp = cal.combine_atomgroups(upper_layers)
        with open(str("_".join(np.unique(ads.resnames)))+'summarystats.csvb', 'w', newline='') as csvfile:
                statWriter = csv.writer(csvfile, delimiter=',')
                for resname in np.unique(ads.resnames): 
                        times, stats = cal.find_adsorption_times(lower_surf_grp,upper_surf_grp,resname,r_c_upper,r_c_lower)
                        statWriter.writerow((resname,np.mean(times),np.std(times)/np.sqrt(len(times))))  
                        hist(times,x_title="Adsorption time\n(number of steps)",y_title="frequency",title=resname)

