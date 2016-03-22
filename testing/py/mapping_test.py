from hybmeshpack import hmscript as hm
import math
global hm, check
hm.check_compatibility("0.3.0")


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


g1 = hm.add_unf_rect_grid([0, 0], [5, 1], 10, 10)
hm.set_boundary_type(g1, 1)
cont = []
for i in range(100):
    x = float(i) / 99.0
    y = 0.1 * math.sin(6 * math.pi * x)
    cont.append([x, y])
for p in cont[::-1]:
    cont.append([p[0], p[1] + 1.0])
cont.append(cont[0])
c1 = hm.create_contour(cont)
hm.set_boundary_type(c1, 2)

print "rectangle to square with sine edges: no, from_contour"
a1 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="no", btypes="from_contour")
check(hm.info_grid(a1) ==
      {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
check(hm.info_contour(a1) ==
      {'btypes': {2: 40}, 'Nnodes': 40, 'subcont': [40], 'Nedges': 40})

print "rectangle to square with sine edges: no, from_grid"
a2 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="no", btypes="from_grid")
check(hm.info_grid(a2) ==
      {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
check(hm.info_contour(a2) ==
      {'btypes': {1: 40}, 'Nnodes': 40, 'subcont': [40], 'Nedges': 40})

print "rectangle to square with sine edges: add_vertices, from_contour"
a3 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="add_vertices", btypes="from_contour")
check(hm.info_grid(a3) ==
      {'cell_types': {11: 4, 4: 80, 13: 4, 14: 4, 15: 8},
       'Nnodes': 313, 'Nedges': 412, 'Ncells': 100})
check(hm.info_contour(a3) ==
      {'btypes': {2: 232}, 'Nnodes': 232, 'subcont': [232], 'Nedges': 232})

print "rectangle to square with sine edges: add_vertices, from_grid"
a4 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="add_vertices", btypes="from_grid")
check(hm.info_grid(a4) ==
      {'cell_types': {11: 4, 4: 80, 13: 4, 14: 4, 15: 8},
       'Nnodes': 313, 'Nedges': 412, 'Ncells': 100})
check(hm.info_contour(a4) ==
      {'btypes': {1: 232}, 'Nnodes': 232, 'subcont': [232], 'Nedges': 232})

print "rectangle to square with sine edges: shift_vertices, from_contour"
a5 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="shift_vertices", btypes="from_contour")
check(hm.info_grid(a5) ==
      {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
check(hm.info_contour(a5) ==
      {'btypes': {2: 40}, 'Nnodes': 40, 'subcont': [40], 'Nedges': 40})

print "rectangle to square with sine edges: shift_vertices, from_grid"
a6 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="shift_vertices", btypes="from_grid")
check(hm.info_grid(a6) ==
      {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
check(hm.info_contour(a6) ==
      {'btypes': {1: 40}, 'Nnodes': 40, 'subcont': [40], 'Nedges': 40})
