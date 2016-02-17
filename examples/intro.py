from HybMeshPyPack import hmscript as hm

# create a base square grid
sqr = hm.AddUnfRectGrid([0, 0], [1, 1], 20, 20)

# create tubes grids
tube_left = hm.AddUnfRectGrid([0, 0], [0.15, 0.5], 10, 20)
tube_right = hm.CopyGeom(tube_left)
tube_bot = hm.CopyGeom(tube_left)
hm.MoveGeom([tube_left, tube_right], 0.425, 1.0)
hm.MoveGeom([tube_bot], 0.425, -0.5)
hm.RotateGeom([tube_left], -45, [0.5, 0.5])
hm.RotateGeom([tube_right], 45, [0.5, 0.5])

# exclude pentagon from base grid
excont = hm.CreateContour([[0.3, 0.6], [0.35, 0.4], [0.65, 0.4],
    [0.7, 0.6], [0.5, 0.7], [0.3, 0.6]])
exsqr = hm.ExcludeContours(sqr, [excont], False)
# ^^^^^^ see fig2

# make grid imposintions with defined buffer sizes
imposed = hm.UniteGrids(exsqr, [(tube_left, 0.05),
    (tube_right, 0.05), (tube_bot, 0.1)])
# ^^^^^^ see fig3

# create boundary grid around pentagon
# since grid domain of `imposed` grid is multiply connected
# we define start_point=end_point which is located on the required
# subcontour of the domain. This excludes all other boundaries from
# the source
bnd1op = hm.BoundaryGridOption(imposed,
        partition=[0, 0.01, 0.02, 0.03],
        direction=1,
        bnd_stepping="keep_shape",
        bnd_step=0.02,
        start_point=[0.3, 0.6],
        end_point=[0.3, 0.6])
bnd1 = hm.BuildBoundaryGrid([bnd1op])

# finally make imposition of the boundary grid
res = hm.UniteGrids(imposed, [(bnd1, 0.05)])
# ^^^^^^ see fig4

# now we can see the result in paraview
hm.ExportVTK(res, "intro.vtk")
