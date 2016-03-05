from hybmeshpack import hmscript

# START OF EXAMPLE

# lower level grid
g1 = hmscript.add_unf_rect_grid([0, 0], [10, 10], 10, 10)
# first imposition grid
g2 = hmscript.add_unf_rect_grid([0, 0], [3, 3], 7, 7)
# second imposition grid
g3 = hmscript.add_unf_circ_grid([5, 5], 1.5, 10, 4)
# impose grids
impgrid = hmscript.unite_grids(g1, [(g2, 2.0), (g3, 1.0)])

# END OF EXAMPLE


print "unite_grids example"
hmscript.export_grid_vtk(impgrid, "_g2.vtk")
if (abs(hmscript.domain_area(impgrid) - 100) > 1e-8):
    raise Exception
if (not hmscript.skewness(impgrid)['ok']):
    raise Exception
