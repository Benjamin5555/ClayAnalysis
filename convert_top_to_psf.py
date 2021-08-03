import parmed as pmd

top = pmd.load_file("control.top")

traj = pmd.load_file("control.trr")

top.save("control.psf",vmd=True)
