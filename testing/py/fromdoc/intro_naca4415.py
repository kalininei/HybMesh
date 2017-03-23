from hybmeshpack import hmscript as hm
hm.check_compatibility("0.5.0", 2)

# === register boundary types
global bfoil, bcirc, binp, bout, bleft, bright, btop, bbot
bfoil = hm.add_boundary_type(1, "foil")
bcirc = hm.add_boundary_type(2, "vortex-generator")
binp = hm.add_boundary_type(3, "input")
bout = hm.add_boundary_type(4, "output")
bleft = hm.add_boundary_type(5, "left")
bright = hm.add_boundary_type(6, "right")
btop = hm.add_boundary_type(7, "top")
bbot = hm.add_boundary_type(8, "bottom")

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
# Here we need a rectangular grid with refinement at foil location.
# First we create two list of doubles defining x and y coordinates,
vert = hm.partition_segment(-1, 1, 0.2, 0.2, [-0.2, 0.03, 0.15, 0.015])
horiz = hm.partition_segment(-1, 3, 0.2, 0.2, [-0.2, 0.03, 1.1, 0.03])
# then build rectangular grid on their basis.
substrate = hm.add_unf_rect_grid(custom_x=horiz, custom_y=vert)


# boundary types for back grid which will be later translated to 3d grid
def _substrate_bfun(x0, y0, x1, y1, bt):
    if abs(x0 - x1) < 1e-12:
        return binp if abs(x0 + 1) < 1e-12 else bout
    else:
        return bbot if abs(y0 + 1) < 1e-12 else btop


hm.set_boundary_type(substrate, bfun=_substrate_bfun)
hm.export_grid_vtk(substrate, "substrate.vtk")

# === build boundary grid around the foil
# First we make custom repartition of the foil to use
# all contour vertices for boundary grid.
# Assign recommended partition distance to points on foil contour:
foilpart = [
    0.01, [0.0, 0.0],     # frontal point
    0.03, [1.0, -0.17],   # backward point
    0.01, [0.25, 0.64],   # upper boundary below vortex generator (finer)
    0.05, [0.5, -0.17],   # lower boundary (coarser)
]
foil = hm.partition_contour(foil, "ref_points", foilpart)

# create boundary grid around foil using contour partition
# obtained above as a horizontal boundary grid partition
foilgrid_bo = hm.BoundaryGridOptions(foil, direction="right",
                                     bnd_stepping="no")
foilgrid_bo.incremental_partition(0.002, 1.5, 7)
foilgrid = hm.build_boundary_grid(foilgrid_bo)


# === boundary grid around vortex generator
vgen = hm.add_circ_contour([0.25, 0.12], 0.015, 32, bcirc)
# build boundary grid
vgen_bo = hm.BoundaryGridOptions(vgen, direction="right",
                                 bnd_step=0.006, bnd_stepping="const")
vgen_bo.incremental_partition(0.002, 1.2, 4)
vgengrid = hm.build_boundary_grid(vgen_bo)
hm.heal_grid(vgengrid)  # get rid of hanging boundary nodes

# === union
res2d = hm.unite_grids(substrate, [(foilgrid, 0.02), (vgengrid, 0.02)], True)

# === wake region
# mesh resolution in the wake region of vortex generator is not fine enough.
# So we apply a patch grid to this region.
vgen_wakegrid = hm.add_unf_rect_grid([0.3, 0.09], [0.4, 0.15], 10, 7)
res2d = hm.unite_grids(res2d, [(vgen_wakegrid, 0.02)], True)

# check if final 2d is fine and stop script if there are bad cells
if not hm.skewness(res2d)['ok']:
    raise Exception("2d grid contains highly skewed cells")

# === 3d grid
# Since we intend to set periodic conditions on z-faces there is no
# need to refine grid towards them. So we use uniform grid in z direction
zmin, zmax, nz = -1.0, 1.0, 20
z_grid = [zmin + i * (zmax - zmin) / nz for i in range(nz + 1)]
res3d = hm.extrude_grid(res2d, z_grid, bleft, bright)

# now we export this grid to fluent msh format keeping periodic conditions
# on left/right and top/bottom surface pairs
hm.export3d_grid_msh(
    res3d, "res3d.msh",
    periodic_pairs=[bleft, bright, [-1, -1, -1], [-1, -1, 1],
                    btop, bbot, [-1, -1, -1], [-1, 3, -1]])
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


def testf():
    from hybmeshpack.hmscript import _dbg as hmdbg
    hmdbg.check_ascii_file(2556162334104726686, "res3d.msh", "dev")


testf()
hm.export3d_grid_vtk(res3d, "g1.vtk")
