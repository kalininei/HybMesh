from hybmeshpack import hmscript as hm
hm.check_compatibility("0.4.0", 2)

# === register boundary types
bfoil = hm.add_boundary_type(1, "foil")
bcirc = hm.add_boundary_type(2, "vortex-generator")
binp = hm.add_boundary_type(3, "input")
bout = hm.add_boundary_type(4, "output")
bleft = hm.add_boundary_type(5, "left")
bright = hm.add_boundary_type(6, "right")

# === build foil contour
# load 2d airfoil coordinates from external file
# to xy float array as [ [x0, y0], [x1, y1], ..... ]
with open('../external_files/naca4415.dat', 'r') as f:
    data = map(float, f.read().split())
    it = iter(data)
    xy = [[x, y] for x, y in zip(it, it)]

# build a contour from read points.
# xy[0] = xy[-1] hence resulting contour will be closed.
# Otherwise we'd have to call xy.append(xy[0][:]) manually.
foil = hm.create_contour(xy, bfoil)

# rotate at 10 degrees
hm.rotate_geom(foil, -10, [0, 0])

# === build substrate grid

# === build boundary grid around the foil
# First we make custom repartition of the foil to use
# all contour vertices for boundary grid.
# Assign recommended partition distance to points on foil contour:
foilpart = [
    0.01, [0.0, 0.0],     # frontal point
    0.03, [1.0, -0.17],   # backward point
    0.01, [0.25, 0.64],   # upper boundary below vortex generator (fine)
    0.05, [0.5, -0.17],   # lower boundary (coarse)
]
foil = hm.partition_contour(foil, "ref_points", foilpart)

# create boundary grid around foil using contour partition
# obtained above as a horizontal boundary grid partition
foilgrid_bo = hm.BoundaryGridOptions(foil, direction="right",
                                     bnd_stepping="no")
foilgrid_bo.incremental_partition(0.002, 1.5, 7)
foilgrid = hm.build_boundary_grid(foilgrid_bo)


# # === boundary grid around vortex generator
# vgen = hm.add_circ_contour([0.25, 0.12], 0.015, 32, bcirc)
# # build boundary grid
# vgen_bo = hm.BoundaryGridOptions(vgen, direction="right",
#                                  bnd_step=0.005, bnd_stepping="const")
# vgen_bo.incremental_partition(0.002, 1.2, 4)
# vgengrid = hm.build_boundary_grid(vgen_bo)

# # === unite boundary grids
# bgrids = hm.unite_grids(foilgrid, [(vgengrid, 0.01)], empty_holes=True)


# hm.export_contour_vtk(foil, "c1.vtk")
# hm.export_contour_vtk(vgen, "c2.vtk")
# hm.export_grid_vtk(foilgrid, "g1.vtk")
# hm.export_grid_vtk(vgengrid, "g2.vtk")
# hm.export_grid_vtk(bgrids, "g3.vtk")
