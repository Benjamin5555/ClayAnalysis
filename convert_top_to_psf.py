import parmed as pmd
import sys

top_path = sys.argv[1]
psf_out = sys.argv[2]

top = pmd.load_file(top_path)


top.save(psf_out,vmd=True)
