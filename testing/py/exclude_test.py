from hybmeshpack import hmscript as hm
global hm, check
hm.check_compatibility("0.2.1")


def check(cond):
    if not cond:
        raise Exception("Test Failed")


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


print "exclude intersecting domain"
g1 = hm.add_unf_circ_grid([0, 0], 1, 20, 10, 1.2, False)
hm.set_boundary_type(g1, 1)
c1 = hm.add_rect_contour([0, 0], [1, 1], [2, 3, 4, 5])
g2 = hm.exclude_contours(g1, c1, True)
g3 = hm.exclude_contours(g1, c1, False)

check_grid(g1, 200, 380, 181, {20: 1, 4: 180})
check_grid(g2, 61, 106, 46, {7: 1, 4: 45})
check_grid(g3, 161, 296, 136, {17: 1, 4: 135})
check_cont(g2, 25, 25, [25], {1: 5, 5: 10, 2: 10})
check_cont(g3, 35, 35, [35], {1: 15, 5: 10, 2: 10})
check_zero(hm.domain_area(g2) + hm.domain_area(g3) - hm.domain_area(g1))

print "exclude internal domain"
g4 = hm.add_unf_rect_grid([0, 0], [0.3, 0.3], 2, 2)
hm.set_boundary_type(g4, 6)
hm.move_geom(g4, -0.5, 0.03)
g5 = hm.exclude_contours(g1, g4, True)
g6 = hm.exclude_contours(g1, g4, False)
[c2] = hm.simplify_contour(g5)
[c3] = hm.simplify_contour(g6)

check_grid(g5, 29, 43, 15, {3: 3, 4: 7, 5: 3, 6: 2})
check_grid(g6, 215, 395, 180, {3: 2, 4: 169, 5: 3, 6: 4, 7: 1, 20: 1})
check_cont(g5, 22, 22, [22], {6: 22})
check_cont(c2, 4, 4, [4], {6: 4})
check_cont(g6, 42, 42, [20, 22], {1: 20, 6: 22})
check_cont(c3, 24, 24, [20, 4], {1: 20, 6: 4})
check_zero(hm.domain_area(g1) - hm.domain_area(c2) - hm.domain_area(c3))

print "exclude external domain"
g7 = hm.add_unf_rect_grid([-10, -10], [10, 10], 1, 1)
g8 = hm.exclude_contours(g1, g7, True)
g9 = hm.exclude_contours(g1, g7, False)

check_zero(hm.domain_area(g1) - hm.domain_area(g8) - hm.domain_area(g9))
check_grid(g8, 200, 380, 181, {20: 1, 4: 180})
check_grid(g9, 0, 0, 0, {})
check_cont(g8, 20, 20, [20], {1: 20})
check_cont(g9, 0, 0, [], {})

print "exclude domain lying within a cell"
c4 = hm.add_rect_contour([-0.02, -0.02], [0.02, 0.02], [6, 7, 8, 9])
g10 = hm.exclude_contours(g1, c4, True)
g11 = hm.exclude_contours(g1, c4, False)

check_zero(hm.domain_area(g1) - hm.domain_area(g10) - hm.domain_area(g11))

check(len(hm.registered_grids()) == 11)
check(len(hm.registered_contours()) == 4)
