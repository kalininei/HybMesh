from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript._dbg import check
global hm, check
hm.check_compatibility("0.4.0")


def check_cont(cont, nn, ne, scont, btypes):
    info = hm.info_contour(cont)
    check(info['Nnodes'] == nn)
    check(info['Nedges'] == ne)
    check(cmp(sorted(info['subcont']), sorted(scont)) == 0)
    for k in btypes.keys():
        check(info['btypes'][k] == btypes[k])

print "contours unite"
g1 = hm.add_unf_rect_grid([-1, -1], [1, 1], 10, 10)
c1 = hm.add_rect_contour([-5, -5], [5, 5], [1, 1, 1, 1])
c2 = hm.unite_contours([c1, g1])
g2 = hm.add_unf_rect_grid([-0.4, -0.4], [0.4, 0.4], 4, 4)
c3 = hm.unite_contours([c2, g2])
c4 = hm.unite_contours([c1, g1, g2])
area_g1 = hm.domain_area(g1)
area_c1 = hm.domain_area(c1)
area_g2 = hm.domain_area(g2)
area_full = hm.domain_area(c4)

check_cont(g1, 40, 40, [40], {0: 40})
check_cont(c1, 4, 4, [4], {1: 4})
check_cont(c2, 44, 44, [4, 40], {0: 40, 1: 4})
check_cont(g2, 16, 16, [16], {0: 16})
check_cont(c3, 60, 60, [4, 16, 40], {0: 56, 1: 4})
check_cont(c4, 60, 60, [4, 16, 40], {0: 56, 1: 4})
check(area_full == 96.64)
check(area_full == area_c1 - area_g1 + area_g2)


print "grid contours to user contours"
c5 = hm.grid_bnd_to_contour(g1, True)
c6 = hm.grid_bnd_to_contour(g2, False)
check_cont(c5, 4, 4, [4], {0: 4})
check_cont(c6, 16, 16, [16], {0: 16})


print "set boundary types"
hm.set_boundary_type(c1, btps=[1, 2, 3, 4])
check_cont(c1, 4, 4, [4], {1: 1, 2: 1, 3: 1, 4: 1})
hm.set_boundary_type(c1, bfun=lambda x0, y0, x1, y1, bo:
                     {1: 0, 2: 0, 3: 0, 4: 15}[bo])
check_cont(c1, 4, 4, [4], {0: 3, 15: 1})
g3 = hm.add_unf_ring_grid([0, 0], 1, 2, 10, 5, 1.4)


def bfun3(x0, y0, x1, y1, bo):
    xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
    if yc > 0:
        return 1
    if (xc * xc + yc * yc) < 1.5:
        return 2
    else:
        return 3


hm.set_boundary_type(g3, bfun=bfun3)
check_cont(g3, 20, 20, [10, 10], {1: 10, 2: 5, 3: 5})

print "separate/simplify"
[c7, c8] = hm.simplify_contour(g3, separate=True)
check_cont(c7, 10, 10, [10], {1: 5, 3: 5})
check_cont(c8, 10, 10, [10], {1: 5, 2: 5})
c9 = hm.create_contour([[0, 0], [3, 1], [6, 0], [9, 0],
                        [9, 3], [8, 6], [9, 9], [9, 12]],
                       [1, 1, 1, 1, 1, 2, 1])
[c10] = hm.simplify_contour(c9, simplify=True, angle=45, separate=True)
check_cont(c10, 5, 4, [4], {1: 3, 2: 1})
c11 = hm.create_contour([[0, 0], [3, 1], [6, 0], [9, 0],
                         [9, 3], [8, 6], [9, 9], [9, 12], [0, 0]],
                        [1, 1, 1, 1, 1, 2, 1])
c12 = hm.create_contour([[-20, -20], [20, -20], [20, 20],
                         [-20, 20], [-20, 0], [-20, -20]], 22)
c13 = hm.unite_contours([c11, c12])
[c14] = hm.simplify_contour(c13, simplify=True, angle=45, separate=False)
[c15, c16] = hm.simplify_contour(c13, simplify=True, angle=45, separate=True)
check_cont(c14, 9, 9, [5, 4], {22: 4, 0: 1, 1: 3, 2: 1})
check_cont(c15, 5, 5, [5], {0: 1, 1: 3, 2: 1})
check_cont(c16, 4, 4, [4], {22: 4})
check(len(hm.registered_contours()) == 16)
hm.remove_all()

print "domain clipping"
c1 = hm.create_contour([[0, 0], [100, 100], [0, 100], [0, 0]], 5)
c2 = hm.add_circ_contour([2, 2], 15, 30, 2)
c3 = hm.clip_domain(c1, c2, "union")
c4 = hm.clip_domain(c1, c2, "difference")
c5 = hm.clip_domain(c1, c2, "intersection")
c6 = hm.clip_domain(c1, c2, "xor")
ar3 = hm.domain_area(c3)
ar4 = hm.domain_area(c4)
ar5 = hm.domain_area(c5)
ar6 = hm.domain_area(c6)
check(abs(ar6 + ar5 - ar3) < 1e-6)
check(abs(ar4 + hm.domain_area(c2) - ar3) < 1e-6)
check_cont(c6, 37, 37, [9, 28], {2: 32, 5: 5})
hm.remove_all()

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10)
c1 = hm.grid_bnd_to_contour(g1, False)
hm.remove_geom(g1)
hm.set_boundary_type(c1, 1)
c2 = hm.add_circ_contour([0.3, 0.3], 0.1, 20, 2)
c3 = hm.add_rect_contour([0.6, 0.5], [0.9, 0.9], 3)
c4 = hm.unite_contours([c1, c2, c3])
check(abs(hm.domain_area(c4) + hm.domain_area(c2) +
          hm.domain_area(c3) - 1) < 1e-6)
c5 = hm.clip_domain(c4, c3, "union", False)
c6 = hm.clip_domain(c5, c2, "union", False)
check(abs(hm.domain_area(c6) - 1.0 < 1e-6))
check(hm.info_contour(c6)['btypes'] == {1: 40})
hm.remove_all()

print "domain reflection"
c1 = hm.create_contour([[0, 0], [1, 0], [2, 1]])
c2 = hm.add_rect_contour([1, 1], [3, 2])
g3 = hm.add_unf_ring_grid([1, 0], 1, 3, 16, 3)
hm.reflect_geom([c1, c2, g3], [0, 0], [1, 1])
hm.remove_all_but([c1])
