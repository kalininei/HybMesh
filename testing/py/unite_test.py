from HybMeshPyPack import hmscript as hm
hm.check_compatibility("0.2.1")


def check(cond):
    import traceback
    if not cond:
        print "TEST FAILED <<<<<<<<<<<<<<<<<<<<<<<"
        traceback.print_stack()


def check_cont(cont, nn, ne, scont, btypes):
    info = hm.info_contour(cont)
    check(info['Nnodes'] == nn)
    check(info['Nedges'] == ne)
    check(cmp(sorted(info['subcont']), sorted(scont)) == 0)
    for k in btypes.keys():
        check(info['btypes'][k] == btypes[k])


def check_grid(grid, nn, ne, nc, ct):
    info = hm.info_grid(grid)
    check(info['Nnodes'] == nn)
    check(info['Nedges'] == ne)
    check(info['Ncells'] == nc)
    for k in ct.keys():
        check(info['cells_types'][k] == ct[k])


def check_zero(a):
    check(abs(a) < 1e-8)

print "unition rect and ring"
g1 = hm.add_unf_rect_grid([-1, -1], [1, 1], 10, 10)
g2 = hm.add_unf_ring_grid([0, 0], 0.5, 0.8, 500, 7)
hm.move_geom([g2], 0, 1)

g3 = hm.unite_grids(g1, [(g2, 0.1)], empty_holes=False)
hm.export_grid_vtk(g3, "g3.vtk")
check(hm.skewness(g3)['ok'])
check(abs(hm.domain_area(g3) - 0.5 * hm.domain_area(g2) - hm.domain_area(g1)) <
      1e-6)


print "two squares: with/without fix_bnd"
g4 = hm.add_unf_rect_grid([10, 10], [15, 15], 5, 5)
g5 = hm.add_unf_rect_grid([10, 10], [11, 11], 30, 30)
hm.move_geom(g5, -0.05, -0.05)
hm.export_grid_vtk(g4, "g4.vtk")
hm.export_grid_vtk(g5, "g5.vtk")
g6 = hm.unite_grids(g5, [(g4, 0.5)], fix_bnd=True)
hm.export_grid_vtk(g6, "g6.vtk")
