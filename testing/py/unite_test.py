from hybmeshpack import hmscript as hm
hm.check_compatibility("0.2.1", 2)


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
        check(info['cell_types'][k] == ct[k])


def check_zero(a):
    check(abs(a) < 1e-8)

print "unition rect and ring"
g1 = hm.add_unf_rect_grid([-1, -1], [1, 1], 10, 10)
g2 = hm.add_unf_ring_grid([0, 0], 0.5, 0.8, 500, 7)
hm.move_geom([g2], 0, 1)

g3 = hm.unite_grids(g1, [(g2, 0.1)], empty_holes=False)
hm.export_grid_vtk(g3, "g3.vtk")
check(hm.skewness(g3)['ok'])
fular = hm.domain_area(g3) - 0.5 * hm.domain_area(g2) - hm.domain_area(g1)
check(abs(fular) < 1e-6)


print "two squares: with/without fix_bnd, bc check"
g4 = hm.add_unf_rect_grid([10, 10], [15, 15], 5, 5)
g5 = hm.add_unf_rect_grid([10, 10], [11, 11], 30, 30)
hm.move_geom(g5, -0.05, -0.05)
hm.set_boundary_type(g4, 3)
hm.set_boundary_type(g5, bfun=lambda x0, y0, x1, y1, b:
                     1 if y0 < 10.8 else 2)
hm.export_grid_vtk(g4, "g4.vtk")
hm.export_grid_vtk(g5, "g5.vtk")
hm.export_contour_vtk(g5, "gc5.vtk")
g6 = hm.unite_grids(g5, [(g4, 0.5)], fix_bnd=False)
check(hm.info_contour(g6)['btypes'] == {1: 5, 2: 1, 3: 18})

g7 = hm.unite_grids(g5, [(g4, 0.5)], fix_bnd=True)
check(hm.info_contour(g7)['btypes'] == {1: 57, 2: 7, 3: 20})


def cosine_boundary(n):
    import math
    ret = []
    for i in range(n):
        x = float(i) / (n - 1)
        ret.append([x, 0.2 * math.sin(2 * math.pi * x)])
    ret.append([0.5, 1])
    ret.append(ret[0])
    return ret

print "unition with a boundary grid"
plist = cosine_boundary(20)
c1 = hm.create_contour(plist, [1] * (len(plist) - 3) + [2, 2])
hm.rotate_geom(c1, 17)
g8 = hm.add_unf_rect_grid([-0.1, -0.1], [1.3, 1.3], 20, 20)
g9 = hm.exclude_contours(g8, c1, "outer")
op1 = hm.BoundaryGridOptions(c1, bnd_stepping="keep_all", bnd_step=0.05)
op1.uniform_partition(0.1, 5)
g10 = hm.build_boundary_grid(op1)
g11 = hm.unite_grids(g9, [(g10, 0.05)])
check(len(hm.skewness(g11)['bad_cells']) <= 3)

print "unition of detached grids"
g12 = hm.add_unf_ring_grid([0, 0], 20, 10, 300, 10)
g13 = hm.add_unf_circ_grid([0, 0], 9.5, 50, 5)
g14 = hm.unite_grids(g12, [(g13, 5)])
check(len(hm.info_contour(g14)['subcont']) == 3)

print "unition of grids with common boundary segments"
g15 = hm.add_unf_circ_grid([0, 0], 1, 16, 5)
g16 = hm.add_unf_rect_grid([-1, -2], [1, -1], 10, 10)
g17 = hm.add_unf_rect_grid([-0.7, -3], [0.3, -2], 20, 20)
g18 = hm.unite_grids(g15, [(g16, 0.3), (g17, 0.3)], fix_bnd=True)
check(hm.skewness(g18)['ok'] and hm.info_grid(g18)['cell_types'][4] == 540)
hm.remove_all()

print "examples for documentation/functionality"
g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10)
g2 = hm.add_unf_rect_grid([0.7, 0.7], [1.3, 1.7], 20, 30)
hm.set_boundary_type(g1, 1)
hm.set_boundary_type(g1, 2)
g3 = hm.unite_grids(g1, [(g2, 0.15)])
g4 = hm.unite_grids(g2, [(g1, 0.15)])
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_contour_vtk(g2, "c2.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
check(hm.skewness(g3)['ok'] and hm.skewness(g4)['ok'])
