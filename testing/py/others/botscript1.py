from hybmeshpack import hmscript as hm

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
hm.export_grid_hmg([cn, mr], "op1.hmg")
cn = hm.unite_grids(cn, [(mr, cn_bufsize)], False, False, 0, "4")

dpath = "others/botscript1_adata/"
files = ["000.msh", "001.msh", "002.msh", "003.msh", "004.msh"]
for f in files:
    print(f)
    dimple = hm.import_grid_gmsh(dpath + f)
    hm.export_grid_hmg([cn, dimple], "op" + f[:3] + ".hmg")
    cn = hm.unite_grids(cn, [(dimple, dm_bufsize)], False, False, 0, "4")

hm.export_grid_vtk(cn, "g1.vtk")
