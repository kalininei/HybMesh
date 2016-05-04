from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript._dbg import check, check_ascii_file
import math
global hm, check
hm.check_compatibility("0.4.2")

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
hm.export_grid_vtk(a1, "g1.vtk")

print "rectangle to square with sine edges: no, from_grid"
a2 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    algo="direct-laplace",
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

print "rectangle grid from 2 straight contours: linear algo"
left_line = hm.create_contour([[0, 0], [0, 1]], 1)
bottom_line = hm.create_contour([[0, 0], [2, 0]], 2)

left_line_part = hm.partition_contour(left_line, "const", 0.2)
bottom_line_part = hm.partition_contour(
    bottom_line, "ref_points",
    [0.2, [0, 0], 0.01, [1, 0], 0.08, [2, 0]])
g1 = hm.add_custom_rect_grid("linear", left_line_part, bottom_line_part)
check(hm.info_contour(g1)['btypes'] == {1: 10, 2: 90})

print "rectangle grid from 4 straight contours: linear algo"
[top_line_part] = hm.copy_geom(bottom_line_part)
hm.move_geom(top_line_part, 0, 1.0)
hm.set_boundary_type(top_line_part, 3)
right_line = hm.create_contour([[2, 0], [2, 1]], 4)
right_line_part = hm.partition_contour(right_line, "ref_points",
                                       [0.1, [2, 0], 0.3, [2, 1]])
g1 = hm.add_custom_rect_grid(
    "linear",
    left_line_part, bottom_line_part,
    right_line_part, top_line_part)
check(hm.info_contour(g1)['btypes'] == {1: 5, 2: 45, 3: 45, 4: 5})

print "linear with disconnected right side"
pts = []
for i in range(100):
    x = 2 + 0.05 * math.sin(4.0 * math.pi * i / 99)
    y = -0.02 + float(i) / 90.0
    pts.append([x, y])
right_line2 = hm.create_contour(pts, 7)
right_line_part2 = hm.partition_contour(right_line2, "const", 0.25)
g1 = hm.add_custom_rect_grid(
    "linear",
    left_line_part, bottom_line_part,
    right_line_part2, top_line_part)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(17821553111849423570, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {1: 5, 2: 45, 3: 45, 7: 5})

print "laplace rectangular grid with linear input"
g1 = hm.add_custom_rect_grid(
    "inverse-laplace",
    left_line_part, bottom_line_part,
    right_line_part, top_line_part)
check(hm.info_contour(g1)['btypes'] == {1: 5, 2: 45, 3: 45, 4: 5})

print "laplace with sined upper"
pts = []
for i in range(100):
    x = 2 * float(i) / 99
    y = 1 + 0.1 * math.sin(8.0 * math.pi * i / 99) - 0.1 * x
    pts.append([x, y])
top_line2 = hm.create_contour(pts, 8)
top_line_part2 = hm.partition_contour(
    top_line2, "ref_points",
    [0.005, [0, 1], 0.23, [2, 1]]
)
left_line_part = hm.partition_contour(left_line, "const", 0.05)
g1 = hm.add_custom_rect_grid(
    "inverse-laplace",
    left_line_part, bottom_line_part,
    None, top_line_part2)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(3760527376004093232, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {8: 45, 1: 40, 2: 45})

g1 = hm.add_custom_rect_grid(
    "orthogonal",
    top_line_part2, left_line_part, bottom_line_part, None)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(8210849007429468971, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {8: 45, 1: 40, 2: 45})

left_line_part2 = hm.partition_contour(
    left_line, "ref_points",
    [0.1, [0, 1], 0.01, [0, 0.5], 0.1, [0, 0]])
g1 = hm.add_custom_rect_grid(
    "direct-laplace",
    top_line_part2, left_line_part2, None, None)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(11365382636530409072, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {8: 90, 1: 52})

hm.export_grid_vtk(g1, "g1.vtk")
hm.export_contour_vtk(g1, "c1.vtk")

print "circ4grid"
g1 = hm.add_circ_rect_grid([10, 10], 10, 1, algo="laplace")
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(11741953394189746279, "g1.vtk")
g1 = hm.add_circ_rect_grid([10, 10], 10, 1, algo="orthogonal-circ")
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(13550186317358967162, "g1.vtk")
