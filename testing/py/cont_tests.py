from HybMeshPyPack import hmscript as hm


def check(cond):
    if not cond:
        raise Exception("Test Failed")


def check_cont(cont, nn, ne, scont, btypes):
    global hm, check
    info = hm.info_contour(cont)
    check(info['Nnodes'] == nn)
    check(info['Nedges'] == ne)
    check(cmp(sorted(info['subcont']), sorted(scont)) == 0)
    for k in btypes.keys():
        check(info['btypes'][k] == btypes[k])

print "contours unite"
g1 = hm.add_unf_rect_grid([-1, -1], [1, 1], 10, 10)
c1 = hm.add_rect_cont([-5, -5], [5, 5], [1, 1, 1, 1])
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
