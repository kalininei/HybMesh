import hybmeshpack.hmscript as hm

outercont = hm.create_contour(
    [[0, 0], [0.2, -0.6], [1, -1.3], [2, -1.3], [2.5, -0.3],
     [2.6, 0.4], [1.5, 1], [0.4, 0.4], [0, 0]])

faultline = hm.create_contour([[0, 1], [1.5, -0.5]])
w1 = [0.64, -0.11]
w2 = [1.54, 0.4]
gw1 = hm.add_unf_circ_grid(w1, 0.05, 16, 5, 1.2, False)
gw2 = hm.add_unf_circ_grid(w2, 0.05, 16, 5, 1.2, False)


faultline = hm.partition_contour(faultline, "const", 0.01, crosses=[outercont])
outercont = hm.matched_partition(outercont, 0.05, 0.5, [faultline])

hm.export_contour_vtk([outercont, faultline], "outer.vtk")

g = hm.pebi_fill(outercont, [faultline])
g = hm.unite_grids(g, [(gw1, 0.01), (gw2, 0.01)], buffer_fill='4')

hm.export_grid_vtk([g], "grid.vtk")
