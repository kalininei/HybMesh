from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript._dbg import check  # NOQA
hm.check_compatibility("0.4.6")

g1 = hm.add_unf_rect_grid(nx=10, ny=10)
z = hm.partition_segment(0, 4, 0.1, 0.1)
g2 = hm.extrude_grid(g1, z, 1, 2, 3)
g3 = hm.tetrahedral_fill(g2)

g4 = hm.add_unf_rect_grid([-1, -1], [2, 2], nx=10, ny=10)
z = hm.partition_segment(-2, 6, 0.4, 0.4)
g5 = hm.extrude_grid(g4, z, 4, 5, 6)

g6 = hm.tetrahedral_fill([g3, g5])

ss = hm.grid3_bnd_to_surface(g6, separate=True)
check(abs(hm.domain_volume(g6) - 68.0) < 1e-8)
check(abs(hm.domain_volume(ss[0]) - 72.0) < 1e-8)
check(abs(hm.domain_volume(ss[1]) - 4.0) < 1e-8)
ss = hm.grid3_bnd_to_surface(g6, separate=False)
check(abs(hm.domain_volume(ss) - 68.0) < 1e-8)

g7 = hm.merge_grids3(g6, g2)
check(abs(hm.domain_volume(g7) - 72.0) < 1e-8)

hm.export3d_grid_vtk(g7, "g7.vtk")
