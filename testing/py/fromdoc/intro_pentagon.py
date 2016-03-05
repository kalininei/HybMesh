from hybmeshpack import hmscript as hm

# create a base square grid
sqr = hm.add_unf_rect_grid([0, 0], [1, 1], 20, 20)

# create tubes grids
tube_left = hm.add_unf_rect_grid([0, 0], [0.15, 0.5], 10, 20)
[tube_right, tube_bot] = hm.copy_geom([tube_left] * 2)
hm.move_geom([tube_left, tube_right], 0.425, 1.0)
hm.move_geom([tube_bot], 0.425, -0.5)
hm.rotate_geom([tube_left], -45, [0.5, 0.5])
hm.rotate_geom([tube_right], 45, [0.5, 0.5])

# exclude pentagon from base grid
excont = hm.create_contour([[0.3, 0.6], [0.35, 0.4], [0.65, 0.4],
                            [0.7, 0.6], [0.5, 0.7], [0.3, 0.6]])
exsqr = hm.exclude_contours(sqr, [excont], "inner")
# ******* see fig2

# make series of grid imposintions with defined buffer sizes
imposed = hm.unite_grids(exsqr, [(tube_left, 0.05),
                                 (tube_right, 0.05),
                                 (tube_bot, 0.1)])
# ******* see fig3

# create boundary grid around pentagon. We use 'imposed' domain as a source
# for boundary grid. Since it is multiply connected
# we define start_point=end_point located on the required
# subcontour of the domain. This excludes all other boundaries from the source.
# Note: here we could have used 'excont' as a source with
#       'direction="right"' option
bnd1op = hm.BoundaryGridOptions(
    imposed,
    partition=[0, 0.01, 0.02, 0.03],
    direction="left",
    bnd_stepping="keep_shape",
    bnd_step=0.02,
    start_point=[0.3, 0.6],
    end_point=[0.3, 0.6])
bnd1 = hm.build_boundary_grid([bnd1op])

# finally make imposition of the boundary grid
res = hm.unite_grids(imposed, [(bnd1, 0.05)])
# ******* see fig4

# now we can see the result in paraview
hm.export_grid_vtk(res, "pentagon.vtk")

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^
if not hm.skewness(res)['ok']:
        raise Exception
