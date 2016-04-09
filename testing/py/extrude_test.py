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
hmdbg.check_ascii_file(5584005863328030835, "s1.vtk")
hmdbg.check_ascii_file(8251898369512273446, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with side bc from contour"
s1 = hm.extrude_grid(g1, [0, 0.1], 0, 0, None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(12347073878777704104, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with non-zero at zfaces and side bc from contour"
s1 = hm.extrude_grid(g1, [-0.1, 0.1, 22], 8, 10, None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(12949446192355444044, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with function at bottom z face"
s1 = hm.extrude_grid(g1, [-0.1, 0.1, 0.2, 0.3],
                     lambda x, y, i: 5 if x < 0.5 else 6,
                     7, None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(11612161261916020638, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with function at bottom and top z face"
s1 = hm.extrude_grid(g1, [0.0, 0.1],
                     lambda x, y, i: 5 if (x < 0.5 and y < 0.6) else 6,
                     lambda x, y, i: i,
                     None)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(4119650101326869963, "s2.vtk")
hm.remove_geom(s1)

print "extrusion with only const boundary types"
z = []
for i in range(0, 101):
    z.append(float(i) / 100)
s1 = hm.extrude_grid(g1, z, 1, 2, 3)
hm.export3d_grid_vtk(s1, None, "s2.vtk")
hmdbg.check_ascii_file(13287298120367059051, "s2.vtk")
hm.remove_geom(s1)
