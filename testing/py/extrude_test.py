from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript import _dbg as hmdbg
hm.check_compatibility("0.4.0")


g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10)
hm.set_boundary_type(g1, bfun=lambda x0, y0, x1, y1, bt: {y0 + y1 == 0: 1,
                                                          x0 + x1 == 2: 2,
                                                          y0 + y1 == 2: 3,
                                                          x0 + x1 == 0: 4}[1])
print "extrusion with zero bc"
s1 = hm.extrude_grid(g1, [0, 0.1, 0.2], 0, 0, 0)
hm.export3d_grid_vtk(s1, "s1.vtk", "s2.vtk")
hmdbg.check_ascii_file(9760293401399903439, "s1.vtk")
hmdbg.check_ascii_file(9133727113521198898, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with side bc from contour"
s1 = hm.extrude_grid(g1, [0, 0.1], 0, 0, None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(17055111601829986000, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with non-zero at zfaces and side bc from contour"
s1 = hm.extrude_grid(g1, [-0.1, 0.1, 22], 8, 10, None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(4721469755399912982, "s2.vtk")

print "extrusion with function at bottom z face"
s1 = hm.extrude_grid(g1, [-0.1, 0.1, 0.2, 0.3],
                     lambda x, y, i: 5 if x < 0.5 else 6,
                     7, None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(5788373827629003185, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with function at bottom and top z face"
s1 = hm.extrude_grid(g1, [0.0, 0.1],
                     lambda x, y, i: 5 if (x < 0.5 and y < 0.6) else 6,
                     lambda x, y, i: i,
                     None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(12326329209190403927, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with only const boundary types"
z = []
for i in range(0, 101):
    z.append(float(i) / 100)
s1 = hm.extrude_grid(g1, z, 1, 2, 3)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(17599784844569501681, "s2.vtk")
hm.remove_geom(s1)

print "full revolution, no contact"
g1 = hm.add_unf_rect_grid([2, 2], [5, 3], 4, 6)
g2 = hm.add_unf_ring_grid([5, 3], 0.1, 0.6, 10, 3)
hm.set_boundary_type(g1, 1)
hm.set_boundary_type(g2, 2)
g3 = hm.unite_grids(g1, [(g2, 0.2)], True)
g4 = hm.revolve_grid(g3, [0, 0], [0.4, 1], 20)
hm.export3d_grid_msh(g4, "g1.msh")
hmdbg.check_ascii_file(11816494025308685888, "g1.msh", "dev")

print "incomplete revolution, with contact, no tri"
g1 = hm.add_unf_circ_grid([4, 4], 15, 32, 15, is_trian=False)
hm.set_boundary_type(g1, 1)
s = 3.5
c1 = hm.add_rect_contour([s, -20], [40, 20], 2)
g1 = hm.exclude_contours(g1, c1, "outer")
g3 = hm.revolve_grid(g1, [s, 0], [s, 1], phi=[0, 10, 20, 30, 40, 50, 60],
                     merge_central=True)
hm.export3d_grid_tecplot(g3, "g1.dat")
hmdbg.check_ascii_file(2839472989048443064, "g1.dat")
hm.export3d_grid_vtk(g3, "g1.vtk")
