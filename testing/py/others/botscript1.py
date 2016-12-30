global hm
from hybmeshpack import hmscript as hm
import os
from os.path import isfile, join

cn_bufsize = 0.100
dm_bufsize = 0.050
mr_min = [-3.500e+0, -5.000e+0]
mr_max = [3.500e+0, 5.000e+0]
mr_step = 0.050

with open("others/botscript1_adata/x1arr.dat") as f:
    x1arr = map(float, f.read().split())
with open("others/botscript1_adata/x2arr.dat") as f:
    x2arr = map(float, f.read().split())

cn = hm.add_unf_rect_grid([0, 0], [1, 1], 3, 3, x1arr, x2arr)
mr = hm.add_unf_rect_grid(mr_min, mr_max, 13, 13, mr_step, mr_step)
cn = hm.unite_grids(cn, [(mr, cn_bufsize)], False, False, 0, "4")

dpath = "others/botscript1_adata/dimples"
files = [join(dpath, f) for f in os.listdir(dpath) if isfile(join(dpath, f))]
print files
for f in files:
    print(f)
    dimple = hm.import_grid_gmsh(f)
    # hm.export_grid_vtk(cn, "a.vtk")
    # hm.export_grid_vtk(dimple, "b.vtk")
    # hm.export_grid_hmg(cn, "a.hmg", fmt="bin")
    # hm.export_grid_hmg(dimple, "b.hmg", fmt="bin")
    cn = hm.unite_grids(cn, [(dimple, dm_bufsize)], False, False, 0, "4")

hm.export_grid_vtk(cn, "g1.vtk")
