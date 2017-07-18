from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript._dbg import check, check_ascii_file, checkdict
global hm, check
hm.check_compatibility("0.4.6")


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
hm.set_boundary_type(
    c1,
    bfun=lambda x0, y0, x1, y1, bo: {1: 0,
                                     2: 0,
                                     3: 0,
                                     4: 15}[bo])
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
[c10] = hm.simplify_contour(c9, simplify=True, angle=60, separate=True)
check_cont(c10, 5, 4, [4], {1: 3, 2: 1})
# hm.export_contour_vtk(c9, "c1.vtk")
# hm.export_contour_vtk(c10, "c2.vtk")
# print hm.info_contour(c10)
# quit()
c11 = hm.create_contour([[0, 0], [3, 1], [6, 0], [9, 0],
                         [9, 3], [8, 6], [9, 9], [9, 12], [0, 0]],
                        [1, 1, 1, 1, 1, 2, 1, 0])
c12 = hm.create_contour([[-20, -20], [20, -20], [20, 20],
                         [-20, 20], [-20, 0], [-20, -20]], 22)
c13 = hm.unite_contours([c11, c12])
[c14] = hm.simplify_contour(c13, simplify=True, angle=60, separate=False)
[c15, c16] = hm.simplify_contour(c13, simplify=True, angle=60, separate=True)
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
hm.remove_all()

print "contour partition"
i1 = hm.create_contour([[0, 0], [1, 0], [2, 1]])
i2 = hm.add_rect_contour([0, 0], [1, 1], [0, 1, 0, 0])
i3 = hm.add_circ_contour([0, 0], 1, 200)
i4 = hm.create_contour([[0, 0], [0.5, 0], [1, 0], [1, 1]])

c1 = hm.partition_contour(i1, "const", 0.18)
c2 = hm.partition_contour(i1, "const", 0.18, 180)
check(hm.info_contour(c1)['Nedges'] == 14)
check(hm.info_contour(c2)['Nedges'] == 13)

c1 = hm.partition_contour(i2, "const", 0.18, 100)
c2 = hm.partition_contour(i2, "const", 0.18, 100, True)
check(hm.info_contour(c1)['Nedges'] == 24)
check(hm.info_contour(c2)['Nedges'] == 23)
hm.export_contour_vtk(c1, "c1.vtk")
hm.export_contour_vtk(c2, "c2.vtk")

c1 = hm.partition_contour(i1, "ref_points", [0.18, [0, 0]])
c2 = hm.partition_contour(i1, "ref_points", [0.01, [0, 0], 0.3, [1, 0]])
c3 = hm.partition_contour(i3, "ref_points", [0.01, [-1, -1],
                                             0.3, [0, 1]])
check(hm.info_contour(c1)['Nedges'] == 14)
check(hm.info_contour(c2)['Nedges'] == 17)
check(hm.info_contour(c3)['Nedges'] == 74)

c1 = hm.partition_contour(i4, "ref_points", [1.0, [0, 0],
                                             1.0, [1, 0],
                                             0.01, [1, 1]])
c2 = hm.partition_contour(i4, "ref_points", [1.0, [0, 0],
                                             1.0, [1, 0],
                                             0.01, [1, 1]], angle0=-1)
check(hm.info_contour(c1)['Nedges'] == 5)
check(hm.info_contour(c2)['Nedges'] == 6)
hm.remove_all()

print "contour partition with given nedges"
cont = hm.create_contour([[0, 5], [1, 5], [0.5, 5.5]])
cont2 = hm.partition_contour(cont, "const", 1, nedges=2)
check(hm.info_contour(cont2)['Nedges'] == 2)
cont2 = hm.partition_contour(cont, "const", 1, nedges=6)
check(hm.info_contour(cont2)['Nedges'] == 6)
hm.export_contour_vtk(cont2, "c.vtk")
check_ascii_file(13813652183113259463, "c.vtk")

ccirc = hm.add_circ_contour([0, 0], 1, 16)
hm.export_contour_vtk(ccirc, "c2.vtk")
cont2 = hm.partition_contour(ccirc, "const", 1, nedges=3, angle0=180.)
hm.export_contour_vtk(cont2, "c3.vtk")
check(hm.info_contour(cont2)['Nedges'] == 3 and
      abs(hm.info_contour(cont2)['length'] - 5.137) < 0.01)

cont2 = hm.partition_contour(ccirc, "const", 1, nedges=22)
check(hm.info_contour(cont2)['Nedges'] == 22)

contcom = hm.unite_contours([ccirc, cont])

cont2 = hm.partition_contour(contcom, "const", 1, nedges=5, angle0=180.)
check(hm.info_contour(cont2)['Nedges'] == 5)

cont2 = hm.partition_contour(contcom, "const", 1, nedges=15)
check(hm.info_contour(cont2)['Nedges'] == 15)

cont2 = hm.partition_contour(
    ccirc, "ref_points", [0.1, [0.5, 0.5], 0.3, [-0.5, -0.5]], nedges=3,
    angle0=180.)
hm.export_contour_vtk(cont2, "c.vtk")
check_ascii_file(3339651189406511573, "c.vtk")

cont2 = hm.partition_contour(
    ccirc, "ref_points", [0.1, [0.5, 0.5], 0.3, [-0.5, -0.5]], nedges=30)
check(hm.info_contour(cont2)['Nedges'] == 30)

cont2 = hm.partition_contour(
    ccirc, "ref_points", [0.1, [0.5, 0.5], 0.3, [-0.5, -0.5]], angle0=0,
    nedges=30)
check(hm.info_contour(cont2)['Nedges'] == 30)
check(hm.domain_area(cont2) == hm.domain_area(ccirc))

hm.export_contour_vtk(ccirc, "c1.vtk")
hm.export_contour_vtk(cont2, "c.vtk")
hm.remove_all()

print "spline tests"
c1 = hm.create_spline_contour([[0, 0], [1, 0.2], [2.4, 0.5], [3, 1]],
                              [0, 5, 0])
check(hm.info_contour(c1)['btypes'] == {0: 56, 5: 44})

c2 = hm.create_spline_contour([[0, 0], [1, 0.2], [2.4, 0.5], [3, 1], [0, 0]],
                              [0, 5, 3, 1])
check(hm.info_contour(c2)['btypes'] == {0: 16, 1: 50, 3: 12, 5: 22})

print "unstructured fill tests"
csqr = hm.add_rect_contour([0, 0], [1, 1], 1)
c1 = hm.partition_contour(csqr, 'ref_points', [0.05, [0, 0], 0.3, [1, 1]])
g1 = hm.triangulate_domain(c1)
checkdict(hm.info_grid(g1), {'cell_types': {3: 124}})

c1 = hm.partition_contour(csqr, 'const', 0.1)
c2 = hm.create_contour([[0.1, 0.2], [0.9, 0.2]])
c2 = hm.partition_contour(c2, 'const', 0.03)
c3 = hm.add_circ_contour([0.3, 0.7], 0.15, 16)
c4 = hm.add_circ_contour([0.7, 0.7], 0.15, 24)
g1 = hm.triangulate_domain(c1, [c2, c3, c4])
check(abs(hm.domain_area(g1)-1)<1e-8);

g1 = hm.triangulate_domain([c1, c3, c4], c2)
check(abs(hm.domain_area(g1) - 0.86123583999)<1e-3);

c1 = hm.create_spline_contour([[-0.2, 0.35], [0.5, 0.1], [1.2, 0.35]])
c1 = hm.partition_contour(c1, "const", 0.05, crosses=[csqr])
c2 = hm.create_contour([[0.3, 0.5], [0.7, 0.5]])
c3 = hm.create_contour([[0.35, 0.4], [0.65, 0.6]])
c4 = hm.create_contour([[0.35, 0.6], [0.65, 0.4]])
c2 = hm.partition_contour(c2, "const", 0.02, crosses=[c3, c4])
c3 = hm.partition_contour(c3, "const", 0.02, crosses=[c2, c4])
c4 = hm.partition_contour(c4, "const", 0.02, crosses=[c2, c3])
c5 = hm.partition_contour(csqr, "const", 0.1, crosses=[c1])
g1 = hm.triangulate_domain(c5, [c1, c2, c3, c4])
# checkdict(hm.info_grid(g1), {'cell_types': {3: 1404}})
check(hm.info_grid(g1)['Ncells'] > 1000 and hm.skewness(g1, 0.75)['ok'])

c1 = hm.partition_contour(csqr, "const", 0.05)
g1 = hm.triangulate_domain(c1, pts=[0.3, [0.3, 0.7],
                                    0.3, [0.7, 0.7],
                                    0.005, [0.5, 0.1]])
# checkdict(hm.info_grid(g1), {'cell_types': {3: 524}})
check(hm.info_grid(g1)['Ncells'] > 500 and hm.skewness(g1)['ok'])

c1 = hm.create_spline_contour([[-0.2, 0.35], [0.5, 0.1], [1.2, 0.35]])
c1 = hm.partition_contour(c1, "const", 0.05, crosses=[csqr])
c2 = hm.partition_contour(csqr, "const", 0.05, crosses=[c1])
g1 = hm.triangulate_domain(c2, c1, pts=[0.3, [0.5, 0.7]])
g2 = hm.triangulate_domain(c2, c1, pts=[0.3, [0.5, 0.7]], fill='4')
# checkdict(hm.info_grid(g1), {'cell_types': {3: 484}})
# checkdict(hm.info_grid(g2), {'cell_types': {3: 4, 4: 225}})
check(hm.skewness(g2)['ok'])
check(hm.skewness(g1)['ok'])

c1 = hm.partition_contour(csqr, "const", 0.05)
c2 = hm.add_circ_contour([0.5, 0.5], 0.3, 120)
g1 = hm.pebi_fill(c1, c2, pts=[0.1, [0.5, 0.5]])
hm.export_grid_vtk(g1, "g1.vtk")
hm.export_contour_vtk([c1, c2], "c1.vtk")
hm.export_contour_vtk(g1, "c2.vtk")
# checkdict(hm.info_grid(g1), {'cell_types': {8: 25, 4: 32, 5: 449,
#                                             6: 747, 7: 233}})
check(hm.skewness(g1)['ok'])


# matched partition
c = hm.add_circ_contour([0, 0], 1, 64)
c1 = hm.create_contour([[-0.65, -0.65], [0.65, 0.65]])
c1 = hm.partition_contour(c1, "const", 0.03)
res = hm.matched_partition(
    c, 0.1, 1.0, [c1], [0.5, [-0.65, 0.65]], angle0=180, power=1)
checkdict(hm.info_contour(res), {'Nedges': 79})

# contour partition
c = hm.add_rect_contour([0, 0], [1, 1], [1, 2, 3, 4])
c1 = hm.partition_contour(
    c, "ref_weights", [0.01, 0, 0.1, 1],
    angle0=180, keep_bnd=True,
    start=[1, 1], end=[0, 0])
hm.export_contour_vtk(c1, "c2.vtk")
check_ascii_file(16637453307356387439, "c2.vtk")

c1 = hm.partition_contour(
    c, "ref_weights", [0.01, 0, 0.1, 1],
    angle0=180, keep_bnd=True,
    start=[1, 1], end=[1, 0])
hm.export_contour_vtk(c1, "c2.vtk")
check_ascii_file(4757620122664042871, "c2.vtk")

c = hm.add_rect_contour([0, 0], [3, 3], [1, 2, 3, 4])
c1 = hm.partition_contour(
    c, "ref_lengths", [0.01, 0.5, 0.1, 0.6, 0.1, 2.4, 0.01, 2.5],
    angle0=180, keep_bnd=True,
    start=[3, 3], end=[0, 3])
hm.export_contour_vtk(c1, "c2.vtk")
check_ascii_file(15116206157781602922, "c2.vtk")

c1 = hm.partition_contour(
    c, "ref_lengths", [0.01, 0, 0.1, 0.5, 0.1, -0.5],
    angle0=180, keep_bnd=True,
    start=[3, 3])
hm.export_contour_vtk(c1, "c2.vtk")
check_ascii_file(10115613561436922700, "c2.vtk")

c1 = hm.partition_contour(
    c, "const", 0.1,
    angle0=180, keep_bnd=True,
    start=[3, 0], end=[3, 3])
hm.export_contour_vtk(c1, "c2.vtk")
check_ascii_file(10535023046019963480, "c2.vtk")


c1 = hm.add_rect_contour([0, 0], [1, 1])
c2 = hm.add_rect_contour([0.5, 0.5], [1.5, 1.5])
c3 = hm.unite_contours([c1, c2])
c4 = hm.decompose_contour(c3)
check(len(c4) == 3)
