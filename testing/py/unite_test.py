from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript._dbg import checkdict
hm.check_compatibility("0.4.5", 2)


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

print "union of rect and ring"
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
check(hm.info_contour(g6)['btypes'] == {1: 3, 2: 1, 3: 18})

# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_grid_vtk(g6, "g6.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")
# hm.export_contour_vtk(g5, "c5.vtk")
# hm.export_contour_vtk(g6, "c6.vtk")
# quit()

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

print "union of with a boundary grid"
plist = cosine_boundary(20)
c1 = hm.create_contour(plist, [1] * (len(plist) - 3) + [2, 2])
hm.rotate_geom(c1, 17)
g8 = hm.add_unf_rect_grid([-0.1, -0.1], [1.3, 1.3], 20, 20)
g9 = hm.exclude_contours(g8, c1, "outer")
op1 = hm.BoundaryGridOptions(c1, bnd_stepping="keep_all", bnd_step=0.05)
op1.uniform_partition(0.1, 5)
g10 = hm.build_boundary_grid(op1)
g11 = hm.unite_grids(g9, [(g10, 0.05)])
check(len(hm.skewness(g11)['bad_cells']) <= 4)

print "union of detached grids"
g12 = hm.add_unf_ring_grid([0, 0], 20, 10, 300, 10)
g13 = hm.add_unf_circ_grid([0, 0], 9.5, 50, 5)
g14 = hm.unite_grids(g12, [(g13, 5)])
check(len(hm.info_contour(g14)['subcont']) == 3)

print "union of grids with common boundary segments"
g15 = hm.add_unf_circ_grid([0, 0], 1, 16, 5)
g16 = hm.add_unf_rect_grid([-1, -2], [1, -1], 10, 10)
g17 = hm.add_unf_rect_grid([-0.7, -3], [0.3, -2], 20, 20)
g18 = hm.unite_grids(g15, [(g16, 0.3), (g17, 0.3)], fix_bnd=True)
check(hm.skewness(g18)['ok'] and hm.info_grid(g18)['cell_types'][4] == 540)

print "large scale difference in non-buffer area"
g19 = hm.add_unf_rect_grid([0, 0], [10, 1], 30, 3)
g20 = hm.add_unf_rect_grid([3, -1], [6, 0], 3, 20)
g21 = hm.unite_grids(g19, [(g20, 0.3)])
check(hm.info_grid(g21)['Ncells'] == 171)
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
hm.remove_all()

g1 = hm.add_unf_circ_grid([0, 1], 0.5, 64, 30, 1.2, False)
hm.set_boundary_type(g1, 1)
g2 = hm.add_unf_rect_grid([0, 0], [1, 1], 30, 30)
hm.set_boundary_type(g2, 2)
g3 = hm.add_unf_rect_grid([-0.3, 0.1], [0.4, 0.7], 30, 40)
hm.set_boundary_type(g3, 3)
g4 = hm.unite_grids(g1, [(g2, 0.1), (g3, 0.1)])
g5 = hm.unite_grids(g1, [(g3, 0.1), (g2, 0.1)])
g6 = hm.unite_grids(g3, [(g2, 0.1), (g1, 0.1)])

# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_contour_vtk(g2, "c2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_grid_vtk(g6, "g6.vtk")

check(hm.skewness(g4)['ok'])
check(hm.skewness(g5)['ok'])
check(len(hm.skewness(g6)['bad_cells']) == 1)
hm.remove_all()

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10)
hm.set_boundary_type(g1, 1)
g21 = hm.add_unf_rect_grid([0.5, 0.5], [1.5, 1.5], 20, 20)
hm.set_boundary_type(g2, 2)
g3 = hm.unite_grids(g1, [(g21, 0.1)])
g22 = hm.add_unf_rect_grid([0.5, 1], [1.5, 2], 20, 20)
hm.set_boundary_type(g22, 2)
g4 = hm.unite_grids(g1, [(g22, 0.1)])
g23 = hm.add_unf_rect_grid([1, 1], [2, 2], 20, 20)
hm.set_boundary_type(g23, 2)
g5 = hm.unite_grids(g1, [(g23, 0.1)])
g24 = hm.add_unf_rect_grid([0.5, 1.01], [1.5, 2.01], 20, 20)
hm.set_boundary_type(g24, 2)
g6 = hm.unite_grids(g1, [(g24, 0.1)])
# hm.move_geom([g3, g4, g5, g6], 3, 0)
# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_grid_vtk(g21, "g21.vtk")
# hm.export_contour_vtk(g21, "c21.vtk")
# hm.export_grid_vtk(g22, "g22.vtk")
# hm.export_contour_vtk(g22, "c22.vtk")
# hm.export_grid_vtk(g23, "g23.vtk")
# hm.export_contour_vtk(g23, "c23.vtk")
# hm.export_grid_vtk(g24, "g24.vtk")
# hm.export_contour_vtk(g24, "c24.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_contour_vtk(g5, "c5.vtk")
# hm.export_grid_vtk(g6, "g6.vtk")
# hm.export_contour_vtk(g6, "c6.vtk")
check(hm.skewness(g3)['ok'])
check(hm.skewness(g4)['ok'])
check(hm.info_grid(g5)['cell_types'].keys() == [4])
check(hm.info_grid(g6)['cell_types'].keys() == [4])
hm.remove_all()


g1 = hm.add_unf_rect_grid([0, 0], [10, 3], 50, 15)
g2 = hm.add_unf_rect_grid([0, -1], [2, 0], 10, 5)
[g3] = hm.copy_geom(g2)
hm.move_geom(g3, 8, 0)
g1 = hm.unite_grids(g1, [(g2, 0), (g3, 0)])
hm.set_boundary_type(g1, 1)
g2 = hm.add_unf_rect_grid([0, -2], [10, -0.5], 15, 3)
hm.set_boundary_type(g2, 2)
hm.rotate_geom([g1, g2], 37)
g3 = hm.unite_grids(g1, [(g2, 0.1)])
g4 = hm.unite_grids(g1, [(g2, 0.5)])
g5 = hm.unite_grids(g1, [(g2, 1)])
g6 = hm.unite_grids(g1, [(g2, 2)])
check(hm.info_grid(g3)['cell_types'][4] == 835)
check(hm.info_grid(g4)['cell_types'][4] == 795)
check(hm.info_grid(g5)['cell_types'][4] == 711)
check(hm.info_grid(g6)['cell_types'][4] == 509)

# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_grid_vtk(g6, "g6.vtk")
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_contour_vtk(g2, "c2.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
hm.remove_all()

g1 = hm.add_unf_rect_grid([0, 0], [4, 4], 20, 20)
hm.set_boundary_type(g1, 1)
g2 = hm.add_unf_rect_grid([3, 3], [5, 5], 10, 10)
hm.set_boundary_type(g2, 2)
g3 = hm.add_unf_circ_grid([2, 2], 1, 6, 3, is_trian=False)
hm.set_boundary_type(g3, 2)
g4 = hm.unite_grids(g1, [(g2, 0)])
g5 = hm.unite_grids(g1, [(g3, 0)])

# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_contour_vtk(g2, "c2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_contour_vtk(g5, "c5.vtk")
check(hm.info_grid(g4)['cell_types'].keys() == [4])
check(hm.info_grid(g5)['cell_types'].keys() == [3, 4, 5, 6, 9, 10])
hm.remove_all()

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10)
g2 = hm.add_unf_rect_grid([0, 0], [0.3, 0.3], 5, 5)
c1 = hm.add_rect_contour([0.06, 0.06], [0.24, 0.24])
g2 = hm.exclude_contours(g2, c1, "inner")
g3 = hm.add_unf_ring_grid([1, 1], 0.3, 0.4, 20, 2)
hm.set_boundary_type(g1, 1)
hm.set_boundary_type(g2, 2)
hm.set_boundary_type(g3, 3)
hm.move_geom(g2, 0.3, 0.3)
g4 = hm.unite_grids(g1, [(g2, 0.1), (g3, 0.1)], empty_holes=False)
g5 = hm.unite_grids(g1, [(g2, 0.1), (g3, 0.1)], empty_holes=True)

# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_contour_vtk(g2, "c2.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")
# hm.export_contour_vtk(g5, "c5.vtk")
check(hm.info_grid(g4)['cell_types'][4] == 113)
check(hm.info_grid(g5)['cell_types'][4] == 110)
check(hm.skewness(g4)['ok'])
check(hm.skewness(g5)['ok'])
hm.remove_all()

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 5, 5)
[g2] = hm.copy_geom([g1])
hm.set_boundary_type(
    g1,
    bfun=lambda x0, y0, x1, y1, bt: 1 if min(x0 + x1, y0 + y1) > 1.2 else 2)
hm.set_boundary_type(
    g2,
    bfun=lambda x0, y0, x1, y1, bt:
    3 if (x0 + x1 == 0 and 0.4 < y0 + y1 < 1.6) else 4)
hm.rotate_geom(g2, 47)
hm.move_geom(g2, 1.16, 0.5)
g3 = hm.unite_grids(g1, [(g2, 0.2)], fix_bnd=True)
g4 = hm.unite_grids(g1, [(g2, 0.2)], fix_bnd=False)

# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_contour_vtk(g2, "c2.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")

check(hm.info_contour(g3)['btypes'] == {1: 2, 2: 16, 3: 2, 4: 17})
check(hm.info_contour(g4)['btypes'] == {2: 17, 4: 17})
hm.remove_all()

g1 = hm.add_unf_ring_grid([0, 0], 3, 6, 256, 15)
c1 = hm.add_rect_contour([-12, 0], [12, 12])
g1 = hm.exclude_contours(g1, c1, "inner")
g2 = hm.add_triangle_grid([-1.4, -5.5], [0, -4.5], [1.4, -5.5], 3)
hm.move_geom(g2, 0, -0.1)
hm.set_boundary_type(g1, 1)
hm.set_boundary_type(g2, 2)
g3 = hm.unite_grids(g1, [(g2, 1)])
g4 = hm.unite_grids(g1, [(g2, 1)], zero_angle_approx=10)
g5 = hm.unite_grids(g1, [(g2, 6.8)], zero_angle_approx=10)
g6 = hm.unite_grids(g1, [(g2, 8)], zero_angle_approx=10)

# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_grid_vtk(g6, "g6.vtk")
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_contour_vtk(g2, "c2.vtk")
# hm.export_contour_vtk(g3, "c3.vtk")
# hm.export_contour_vtk(g4, "c4.vtk")

n1 = hm.info_contour(g3)['Nedges']
n2 = hm.info_contour(g4)['Nedges']
n3 = hm.info_contour(g5)['Nedges']
n4 = hm.info_contour(g6)['Nedges']
check(n1 == 286 and n2 < n1 and n3 < n2 and n4 < n3)
check(hm.skewness(g6)['ok'])

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 20, 20)
g2 = hm.add_unf_rect_grid([0.3, 0.3], [0.6, 0.6], 7, 7)
g3 = hm.unite_grids(g1, [(g2, 0.1)], buffer_fill="3")
g4 = hm.unite_grids(g1, [(g2, 0.1)], buffer_fill="4")
checkdict(hm.info_grid(g3), {'cell_types': {3: 196, 4: 349}})
# checkdict(hm.info_grid(g4), {'cell_types': {4: 435}})
check(hm.skewness(g4)['ok'])
