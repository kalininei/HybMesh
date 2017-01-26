from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript._dbg import check, check_ascii_file, checkdict
import math
global hm, check, checkdict, math
hm.check_compatibility("0.4.3")

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
checkdict(
    hm.info_grid(a1),
    {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
checkdict(
    hm.info_contour(a1),
    {'btypes': {2: 40}, 'Nnodes': 40, 'subcont': [40], 'Nedges': 40})
hm.export_grid_vtk(a1, "g1.vtk")

print "rectangle to square with sine edges: no, from_grid"
a2 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    algo="direct_laplace",
    snap="no", btypes="from_grid")
checkdict(
    hm.info_grid(a2),
    {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
checkdict(
    hm.info_contour(a2),
    {'btypes': {1: 40}, 'Nnodes': 40, 'subcont': [40], 'Nedges': 40})


print "rectangle to square with sine edges: add_vertices, from_contour"
a3 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="add_vertices", btypes="from_contour")
checkdict(
    hm.info_grid(a3),
    {'cell_types': {11: 4, 4: 80, 13: 4, 14: 4, 15: 8},
     'Nnodes': 313, 'Nedges': 412, 'Ncells': 100})
checkdict(
    hm.info_contour(a3),
    {'btypes': {2: 232}, 'Nnodes': 232, 'subcont': [232], 'Nedges': 232})

print "rectangle to square with sine edges: add_vertices, from_grid"
a4 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="add_vertices", btypes="from_grid")
checkdict(
    hm.info_grid(a4),
    {'cell_types': {11: 4, 4: 80, 13: 4, 14: 4, 15: 8},
     'Nnodes': 313, 'Nedges': 412, 'Ncells': 100})
checkdict(
    hm.info_contour(a4),
    {'btypes': {1: 232}, 'Nnodes': 232, 'subcont': [232], 'Nedges': 232})

print "rectangle to square with sine edges: shift_vertices, from_contour"
a5 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="shift_vertices", btypes="from_contour")
checkdict(
    hm.info_grid(a5),
    {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
checkdict(
    hm.info_contour(a5),
    {'btypes': {2: 40}, 'Nnodes': 40, 'subcont': [40], 'Nedges': 40})

print "rectangle to square with sine edges: shift_vertices, from_grid"
a6 = hm.map_grid(
    g1, c1,
    [[0, 0], [5, 0], [5, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    snap="shift_vertices", btypes="from_grid")
checkdict(
    hm.info_grid(a6),
    {'cell_types': {4: 100}, 'Nnodes': 121, 'Nedges': 220, 'Ncells': 100})
checkdict(
    hm.info_contour(a6),
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
right_line_part2 = hm.partition_contour(
    right_line2, "const", 0.25, angle0=180.,
    nedges=hm.info_contour(left_line_part)['Nedges'])
g1 = hm.add_custom_rect_grid(
    "linear",
    left_line_part, bottom_line_part,
    right_line_part2, top_line_part)
hm.export_grid_vtk(g1, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {1: 5, 2: 45, 3: 45, 7: 5})

print "laplace rectangular grid with linear input"
g1 = hm.add_custom_rect_grid(
    "inverse_laplace",
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
    [0.005, [0, 1], 0.23, [2, 1]],
    angle0=180.
)
left_line_part = hm.partition_contour(left_line, "const", 0.05)
g1 = hm.add_custom_rect_grid(
    "inverse_laplace",
    left_line_part, bottom_line_part,
    None, top_line_part2)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(16985162519033738007, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {8: 45, 1: 40, 2: 45})

g1 = hm.add_custom_rect_grid(
    "orthogonal",
    top_line_part2, left_line_part, bottom_line_part, None)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(17680407579251078528, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {8: 45, 1: 40, 2: 45})

left_line_part2 = hm.partition_contour(
    left_line, "ref_points",
    [0.1, [0, 1], 0.01, [0, 0.5], 0.1, [0, 0]])
g1 = hm.add_custom_rect_grid(
    "direct_laplace",
    top_line_part2, left_line_part2, None, None)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(8646985382005859307, "g1.vtk")
check(hm.info_contour(g1)['btypes'] == {8: 90, 1: 52})


print "circ4grid"
g1 = hm.add_circ_rect_grid([10, 10], 10, 1, algo="laplace")
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(15941694762539007875, "g1.vtk")
g1 = hm.add_circ_rect_grid([10, 10], 10, 1, algo="orthogonal_circ")
hm.export_grid_vtk(g1, "g2.vtk")
check_ascii_file(16442454580617227043, "g2.vtk")

print "examples for 'Map grid' documentation"
hm.remove_all()
g1 = hm.add_unf_rect_grid([0, 0], [2, 1], 30, 20)
c1 = hm.add_rect_contour([0, 0], [2, 1], [1, 2, 3, 4])
g1 = hm.exclude_contours(g1, c1, "outer")


def addline(c, p1, p2, a, k):
    vec = [p2[0] - p1[0], p2[1] - p1[1]]
    vlen = math.sqrt(vec[0]**2 + vec[1]**2)
    vec[0] /= vlen
    vec[1] /= vlen
    for i in range(100):
        t = float(i) / 100
        m = a * math.sin(2.0 * math.pi * k * t)
        p = [t * s2 + (1 - t) * s1 for (s1, s2) in zip(p1, p2)]
        perpp = [-vec[1] * m, vec[0] * m]
        p = [p[0] + perpp[0], p[1] + perpp[1]]
        c.append(p)

p1 = [3, 0]
p2 = [4, -0.3]
p3 = [4.2, 0.3]
p4 = [2.9, 0.8]
cont1 = []
addline(cont1, p1, p2, 0.02, 3.0)
addline(cont1, p2, p3, 0.02, 2.0)
addline(cont1, p3, p4, 0.05, 1.0)
addline(cont1, p4, p1, 0.02, 3.0)
cont1.append(cont1[0])
c = hm.create_contour(cont1)
g2 = hm.map_grid(g1, c, [[0, 0], [2, 0], [2, 1], [0, 1]], [p1, p2, p3, p4],
                 algo="direct_laplace")
g3 = hm.map_grid(
    g1, c, [[0, 0], [2, 0], [2, 1], [0, 1], [1, 0]],
    [p1, p2, p3, p4, [3.14, 0]],
    algo="direct_laplace")

g4 = hm.map_grid(
    g1, c, [[0, 0], [2, 0], [2, 1], [0, 1], [0, 0.5], [2, 0.5]],
    [p1, p2, p3, p4, [2.9, 0.66], [4.12, 0.12]],
    algo="direct_laplace")

g5 = hm.map_grid(
    g1, c, [[0, 0]], [p1], algo="direct_laplace")

# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_contour_vtk(c, "ctarget.vtk")
# hm.export_contour_vtk(g1, "cbase.vtk")
# hm.export_grid_vtk(g5, "gtar.vtk")
# hm.export_contour_vtk(g5, "ctar.vtk")
hm.remove_all()

g1 = hm.add_unf_ring_grid([0, 0], 1, 1.5, 32, 12, 0.95)
c1 = hm.add_circ_contour([0, 0], 1, 100)
cm = hm.add_rect_contour([-30, -30], [0, 30])
c1 = hm.clip_domain(c1, cm, "intersection")
hm.scale_geom(c1, 300, 100)
cm = hm.add_circ_contour([0, 0], 1, 100)
c1 = hm.clip_domain(c1, cm, "union")
c2 = hm.add_circ_contour([-2, 0], 0.3, 100)
hm.scale_geom(c2, 100, 30)
c = hm.unite_contours([c1, c2])

gtmp = hm.map_grid(g1, c, [[-1.5, 0], [-1, 0]], [[-3, 0], [-2, 0.01]],
                   algo="direct_laplace", return_invalid=True)
hm.export_grid_vtk(gtmp, "gtmp.vtk")

g2 = hm.map_grid(g1, c, [[-1.5, 0], [-1, 0]], [[-3, 0], [-2.3, 0]],
                 algo="direct_laplace")


g3 = hm.map_grid(g1, c, [[-1.5, 0], [-1, 0], [0, 1.5], [0, -1.5]],
                 [[-3, 0], [-2.3, 0.0], [-2, 0.73], [-2, -0.73]],
                 algo="direct_laplace")

# hm.export_grid_vtk(g1, "g.vtk")
# hm.export_grid_vtk(g2, "g1.vtk")
# hm.export_grid_vtk(g3, "g2.vtk")
hm.remove_all()

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10)
c1 = hm.create_contour([[0, 0], [1, 0], [1, 1], [0.5, 0.2], [0, 1], [0, 0]])
g2 = hm.map_grid(
    g1, c1,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    algo="direct_laplace", return_invalid=True)
g3 = hm.map_grid(
    g1, c1,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    algo="inverse_laplace", return_invalid=True)
g4 = hm.map_grid(
    g3, g1,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    algo="direct_laplace", return_invalid=True)
g5 = hm.map_grid(
    g3, g1,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    algo="inverse_laplace", return_invalid=True)
# hm.export_contour_vtk(g1, "c1.vtk")
# hm.export_contour_vtk(g3, "c2.vtk")
# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
hm.remove_all()

cont = []
for i in range(0, 100):
    t = float(i) / 99
    y = -t if t < 0.5 else t - 1
    cont.append([10 * t - 5, 10 * y])
cont.extend([[5, 4], [-5, 4], cont[0]])
c = hm.create_contour(cont)

g = hm.add_unf_rect_grid([0, 0], [1, 1], 40, 20)
g1 = hm.map_grid(
    g, c,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [[-5, 0], [5, 0], [5, 4], [-5, 4]],
    algo="direct_laplace")
g2 = hm.map_grid(
    g, c,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [[-5, 0], [5, 0], [5, 4], [-5, 4]],
    algo="inverse_laplace")
# hm.export_grid_vtk(g, "g.vtk")
# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_contour_vtk(g, "c1.vtk")
# hm.export_contour_vtk(g1, "c2.vtk")
hm.remove_all()

g = hm.add_unf_rect_grid([0, 0], [2, 2], 40, 40)
c = hm.add_rect_contour([0.5, 0.5], [1.5, 1.5])
g = hm.exclude_contours(g, c, "inner")

c = hm.add_circ_contour3([0, 2], [0, 0], 0.7, 64)
c2 = hm.add_rect_contour([-100, -100], [0, 100])
c = hm.clip_domain(c, c2, "intersection")
c2 = hm.add_rect_contour([0, 0], [2, 2])
c = hm.clip_domain(c, c2, "union")
c2 = hm.add_rect_contour([0.2, 0.9], [1.5, 1.4])
c3 = hm.add_circ_contour3([1.5, 0.9], [0.2, 0.9], 1, 64)
c2 = hm.clip_domain(c2, c3, "difference")
c = hm.clip_domain(c, c2, "difference")
g1 = hm.map_grid(
    g, c,
    [[0, 0], [0, 2], [0.5, 0.5], [1.5, 0.5], [1.5, 1.5], [0.5, 1.5]],
    [[0, 0], [0, 2], [0.2, 0.9], [1.5, 0.9], [1.5, 1.4], [0.2, 1.4]],
    algo="inverse_laplace")
g2 = hm.map_grid(
    g, c,
    [[0, 0], [0, 2], [0.5, 0.5], [1.5, 0.5], [1.5, 1.5], [0.5, 1.5]],
    [[0, 0], [0, 2], [0.2, 0.9], [1.5, 0.9], [1.5, 1.4], [0.2, 1.4]],
    algo="direct_laplace")
# hm.export_contour_vtk(c, "c.vtk")
# hm.export_contour_vtk(g, "c1.vtk")
# hm.export_grid_vtk(g, "g.vtk")
# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
hm.remove_all()

cont = []
addline(cont, [1, 1], [0, 1], 0.35, 1)
cont.extend([[0, 1], [0, 0], [1, 0], [1, 1]])
c1 = hm.create_contour(cont)
c1 = hm.partition_contour(c1, "const", 0.08, angle0=10)

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 7, 9)
g2 = hm.map_grid(g1, c1, [[0, 1], [1, 1]], [[0, 1], [1, 1]])
g3 = hm.map_grid(g1, c1, [[0, 1], [1, 1]], [[0, 1], [1, 1]],
                 snap="shift_vertices")
g4 = hm.map_grid(g1, c1, [[0, 1], [1, 1]], [[0, 1], [1, 1]],
                 snap="add_vertices")

# hm.export_contour_vtk(c1, "c1.vtk")
# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
hm.remove_all()

print "examples for custom rectangular grid documentation"


def addline(c, p1, p2, a, k):
    vec = [p2[0] - p1[0], p2[1] - p1[1]]
    vlen = math.sqrt(vec[0]**2 + vec[1]**2)
    vec[0] /= vlen
    vec[1] /= vlen
    for i in range(100):
        t = float(i) / 99
        m = a * math.sin(2.0 * math.pi * k * t)
        p = [t * s2 + (1 - t) * s1 for (s1, s2) in zip(p1, p2)]
        perpp = [-vec[1] * m, vec[0] * m]
        p = [p[0] + perpp[0], p[1] + perpp[1]]
        c.append(p)

pnt = []
addline(pnt, [0, 0], [0.1, 1], 0.05, 0.5)
leftc = hm.create_contour(pnt, 1)
pnt = []
addline(pnt, [0, 0], [2, 0], 0.02, 1)
botc = hm.create_contour(pnt, 2)
pnt = []
addline(pnt, [2, 0], [1.93, 1.03], 0.05, 0.5)
rightc = hm.create_contour(pnt, 3)
pnt = []
addline(pnt, [0.1, 1], [1.93, 1.03], 0.06, 1)
topc = hm.create_contour(pnt, 4)
fullc = hm.unite_contours([leftc, botc, rightc, topc])

left = hm.partition_contour(
    leftc, "ref_points",
    [0.2, [0, 0], 0.2, [0.1, 1], 0.05, [0, 0.5]])
right = hm.partition_contour(rightc, "const", 0.12)
top = hm.partition_contour(topc, "const", 0.1)
bot = hm.partition_contour(botc, "const", 0.11)

fullc1 = hm.unite_contours([left, bot, right, top])
g1 = hm.add_custom_rect_grid("linear", left, bot, right, top)


pnt = []
addline(pnt, [0, 0], [0.1, 1], 0.1, 0.5)
leftc = hm.create_contour(pnt, 1)
pnt = []
addline(pnt, [2, 0], [1.93, 1.03], 0.05, 1.5)
rightc = hm.create_contour(pnt, 3)
pnt = []
addline(pnt, [0, 0], [2, 0], 0.08, 1)
botc = hm.create_contour(pnt, 2)
pnt = []
addline(pnt, [0.1, 1], [1.93, 1.03], -0.08, 1)
topc = hm.create_contour(pnt, 4)
fullc = hm.unite_contours([leftc, botc, rightc, topc])

left = hm.partition_contour(
    leftc, "ref_points",
    [0.2, [0, 0], 0.2, [0.1, 1], 0.05, [0, 0.5]],
    angle0=180.)
right = hm.partition_contour(rightc, "const", 0.12, angle0=180.)
top = hm.partition_contour(topc, "const", 0.095, angle0=180.)
bot = hm.partition_contour(
    botc, "ref_points",
    [0.03, [0, 0], 0.15, [2, 0], 0.15, [1, 0]], angle0=180.)

fullc2 = hm.unite_contours([left, bot, right, top])
g2 = hm.add_custom_rect_grid("direct_laplace", left, bot, right, top)


hm.export_contour_vtk(fullc1, "left1.vtk")
hm.export_grid_vtk(g1, "g1.vtk")
hm.export_contour_vtk(fullc2, "left2.vtk")
hm.export_grid_vtk(g2, "g2.vtk")

pnt = []
addline(pnt, [0, 0], [0.1, 1], 0.1, 0.5)
leftc = hm.create_contour(pnt, 1)
pnt = []
addline(pnt, [2, 0], [1.93, 1.03], 0.05, 1.5)
rightc = hm.create_contour(pnt, 3)
pnt = []
addline(pnt, [0, 0], [2, 0], 0.08, 1)
botc = hm.create_contour(pnt, 2)
pnt = []
addline(pnt, [0.1, 1], [1.93, 1.03], -0.18, 1)
topc = hm.create_contour(pnt, 4)
fullc = hm.unite_contours([leftc, botc, rightc, topc])

left = hm.partition_contour(
    leftc, "ref_points",
    [0.2, [0, 0], 0.2, [0.1, 1], 0.05, [0, 0.5]],
    angle0=180.)
right = hm.partition_contour(rightc, "const", 0.12, angle0=180.)
top = hm.partition_contour(topc, "const", 0.095, angle0=180.)
bot = hm.partition_contour(
    botc, "ref_points",
    [0.03, [0, 0], 0.14, [2, 0], 0.15, [1, 0]], angle0=180.)

fullc3 = hm.unite_contours([left, bot, right, top])
g3 = hm.add_custom_rect_grid("inverse_laplace", left, bot, right, top)

pnt = []
addline(pnt, [0, 0], [0.1, 1], 0.1, 0.5)
leftc = hm.create_contour(pnt, 1)
pnt = []
addline(pnt, [2, 0], [1.93, 1.03], 0.05, 1.5)
rightc = hm.create_contour(pnt, 3)
pnt = []
addline(pnt, [0, 0], [2, 0], 0.08, 1)
botc = hm.create_contour(pnt, 2)
pnt = []
addline(pnt, [0.1, 1], [1.93, 1.03], -0.18, 1)
topc = hm.create_contour(pnt, 4)
fullc = hm.unite_contours([leftc, botc, rightc, topc])

left = hm.partition_contour(
    leftc, "ref_points",
    [0.2, [0, 0], 0.2, [0.1, 1], 0.05, [0, 0.5]])
right = hm.partition_contour(rightc, "const", 0.12)
top = hm.partition_contour(topc, "const", 0.095)
bot = hm.partition_contour(
    botc, "ref_points",
    [0.03, [0, 0], 0.14, [2, 0], 0.15, [1, 0]])

fullc4 = hm.unite_contours([left, bot, right, top])
g4 = hm.add_custom_rect_grid("orthogonal", left, bot, right, top)

left = hm.partition_contour(
    leftc, "ref_points",
    [0.05, [0, 0], 0.05, [0.1, 1], 0.2, [0, 0.5]],
    angle0=180.)
right = hm.partition_contour(rightc, "const", 0.12, angle0=180.)
top = hm.partition_contour(topc, "const", 0.095, angle0=180.)
bot = hm.partition_contour(
    botc, "ref_points",
    [0.03, [0, 0], 0.14, [2, 0], 0.15, [1, 0]],
    angle0=180.)

fullc5 = hm.unite_contours([left, bot, right, top])
g5 = hm.add_custom_rect_grid("linear_tfi", left, bot, right, top)
g6 = hm.add_custom_rect_grid("hermite_tfi", left, bot, right, top,
                             [0.4] * 4)

# hm.export_contour_vtk(fullc1, "left1.vtk")
# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_contour_vtk(fullc2, "left2.vtk")
# hm.export_grid_vtk(g2, "g2.vtk")
# hm.export_contour_vtk(fullc3, "left3.vtk")
# hm.export_grid_vtk(g3, "g3.vtk")
# hm.export_contour_vtk(fullc4, "left4.vtk")
# hm.export_grid_vtk(g4, "g4.vtk")
# hm.export_contour_vtk(fullc5, "left1.vtk")
# hm.export_grid_vtk(g5, "g5.vtk")
# hm.export_grid_vtk(g6, "g6.vtk")
hm.remove_all()

left = hm.create_contour([[0, 0], [0, 1]])
right = hm.create_contour([[1, 0], [1, 1]])
top = hm.create_contour([[0, 1], [1, 1]])
pnt = []
addline(pnt, [0, 0], [1, 0], 0.16, 0.5)
bot = hm.create_contour(pnt, 2)
left = hm.partition_contour(
    left, "ref_points",
    [0.003, [0, 0], 0.1, [0, 1]], nedges=20)
right = hm.partition_contour(
    right, "ref_points",
    [0.003, [1, 0], 0.1, [1, 1]], nedges=20)
top = hm.partition_contour(top, "ref_points", [1, [0, 1], 1, [1, 1]],
                           nedges=10)
bot = hm.partition_contour(bot, "ref_points", [1, [0, 0], 1, [1, 0]],
                           nedges=10)

fullc = hm.unite_contours([left, bot, right, top])
g1 = hm.add_custom_rect_grid("hermite_tfi", left, bot, right, top,
                             [0.0, 3, 0.0, 0.0])

# hm.export_grid_vtk(g1, "g1.vtk")
# hm.export_contour_vtk(fullc, "c1.vtk")
hm.remove_all()
