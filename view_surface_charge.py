from clayAnalysis import ClayAnalysis 
import numpy as np
import MDAnalysis as mda

u = mda.Universe("../Caesium/topol.tpr","../Caesium/clay_ions_eq_d.gro")

cal = ClayAnalysis(u)

surfaces = cal.generate_surface_group("waters")
upper_surf = surfaces[0]

ag_upper = upper_surf[0]
for i in upper_surf[1:]:
    ag_upper = ag_upper + i

surfaces = cal.generate_surface_group("mineral")
upper_surf = surfaces[0]
for i in upper_surf[1:]:
    ag_upper = ag_upper + i

cal.plot_surface_charge(ag_upper,2)
