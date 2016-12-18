from hybmeshpack import hmscript as hm
import copy  # NOQA
import math  # NOQA
global hm, check  # NOQA
hm.check_compatibility("0.4.6")

# options
step_main = 20.
step_bl = 2.
bstep = [step_bl, 0, step_main, 0.1, step_main, 0.9, step_bl, 1]
bstep1 = [step_bl, 0, step_main, 0.1]
bstep_main = [step_bl, 0, step_main, 0.05, step_main, 0.95, step_bl, 1]

# importing
cont = hm.import_contour_hmc("../../external_files/RRJ2.hmc")
hm.export_contour_vtk(cont, "cont.vtk")

# =characteristic points
# lower bnd
c = [None] * 100
c[0] = [-906.994, -42.374]
c[1] = [2099.83, -33.8302]
c[2] = [1913.45, 1508.89]
c[3] = [-729.668, 1501.38]
# upper bnd
c[5] = [271.595, 1587.21]
c[6] = [1124.28, 1589.55]
c[7] = [1237.86, 2013.14]
c[8] = [156.175, 2009.8]
# left appendix line
c[9] = [271.027, 1592.65]
c[10] = [292.402, 1646.47]
# right appendix line
c[11] = [1124.81, 1594.99]
c[12] = [1103.15, 1649.06]
# left top zone rb point
c[13] = [201.892, 1966.72]
# right top zone lb point
c[14] = [1190.18, 1969.44]
# left input
c[15] = [-588.749, 1472.73]
c[16] = [-581.982, 1472.74]
# right input
c[17] = [1767.55, 1479.88]
c[18] = [1773.05, 1479.43]
# left top zone
c[19] = [161.65, 1961.3]
c[20] = [156.118, 1961.29]
c[21] = [110.831, 1963.87]
c[22] = [112.461, 1968.8]
c[23] = [129.461, 1990.69]
c[32] = [134.915, 1995.5]
# left chairs
c[24] = [-535.8674316, 1081.113525]
c[25] = [-338.0876465, 1080.345825]
c[26] = [-21.50970459, 1080.312256]
c[27] = [172.5718689, 1079.886475]
c[28] = [-459.754, 593.154]
c[29] = [-410.738, 593.289]
c[30] = [52.5625, 594.558]
c[31] = [101.578, 594.693]
# right top zone
c[33] = [1234.68, 1956.28]
c[34] = [1268.93, 1990]
c[35] = [1280.92, 1972.36]
c[36] = [1283.32, 1967.14]
c[37] = [1235.21, 1926.42]
c[38] = [1265.23, 1926.5]
c[39] = [1281.34, 1955.58]
# right chairs
c[40] = [1554.99, 1081.6]
c[41] = [1743.26, 1080.36]
c[42] = [1626.44, 598.872]
c[43] = [1675.87, 599.006]
c[44] = [2069.82, 1079.89]
c[45] = [2138.76, 600.276]
c[46] = [2187.77, 600.41]
c[47] = [2147.26, 457.037]
c[48] = [2188.22, 457.037]
c[49] = [2147.59, 302.028]

for i, p in enumerate(c):
    if p is not None:
        c[i] = hm.get_point(cont, vclosest=p)

# closed contours assembling
cn = [cont]
cn.append(hm.create_contour([c[9], c[10]]))
cn.append(hm.create_contour([c[11], c[12]]))
cn.append(hm.create_contour([c[11], c[12]]))
cn.append(hm.create_contour([c[8], c[13]]))
cn.append(hm.create_contour([c[7], c[14]]))
cn.append(hm.create_contour([c[13], c[14]]))
cn.append(hm.create_contour([c[5], c[6]]))
d = hm.unite_contours(cn)
dcm = hm.decompose_contour(d)
c_ltop = hm.pick_contour(c[21], dcm)
c_rtop = hm.pick_contour(c[38], dcm)
c_lhole = hm.pick_contour([140, 1936], dcm)
c_rhole = hm.pick_contour([1254, 1940], dcm)
c_top = hm.pick_contour([687, 2096], dcm)
c_mid = hm.pick_contour([282, 1803], dcm)
c_lappendix = hm.pick_contour([252, 1617], dcm)
c_rappendix = hm.pick_contour([1164, 1611], dcm)
c_bot = hm.pick_contour(c[0], dcm)
c_lchairs = hm.pick_contour(c[24], dcm)
c_rchairs = hm.pick_contour([1390, 1105], dcm)

# bottom grid: g1
[bc, rc, t1, t2, t3, t4, t5, t6, t7, lc] = hm.extract_subcontours(
    c_bot, [c[i] for i in [0, 1, 2, 18, 17, 6, 5, 16, 15, 3, 0]])
t2 = hm.create_contour([c[18], c[17]])
t6 = hm.create_contour([c[16], c[15]])
tc = hm.unite_contours([t1, t2, t3, t4, t5, t6, t7])

lc = hm.partition_contour(lc, "ref_weights", bstep_main, 90, start=c[0])
tc = hm.partition_contour(tc, "ref_weights", bstep_main, 90, start=c[2],
                          keep_pts=[c[5], c[6]])
g1 = hm.add_custom_rect_grid("orthogonal", lc, tc, rc, bc)
g1cont = hm.grid_bnd_to_contour(g1, False)

# left appendix: g4
bs = [0, step_bl, 2 * step_bl, 3 * step_bl]
bgo = hm.BoundaryGridOptions(c_lappendix, bs, "left", 3 * step_bl,
                             "const", start_point=c[10], end_point=c[9])
g41 = hm.build_boundary_grid(bgo)
lc = hm.create_contour([hm.get_point(g41, 104), hm.get_point(g41, 39)])
lc = hm.partition_contour(lc, "const", 4 * step_bl)
u = hm.unite_contours([g41, lc])
d = hm.decompose_contour(u)
lc = hm.pick_contour([281, 1612.484], d)
g4 = hm.triangulate_domain(lc, fill='4')
g4 = hm.unite_grids(g41, [(g4, 0)])
hm.heal_grid(g4, simplify_boundary=50)

# right appendix: g5
g5 = hm.map_grid(g4, c_rappendix,
                 [c[9], c[10], [231.5, 1608.5], [286.399, 1647.06]],
                 [c[11], c[12], [1164, 1611], [1109.88, 1649.44]],
                 "no", "vertex", is_reversed=True)

# left input: g6
p1 = hm.get_point(c_bot, vclosest=[-592.302, 1483.16])
p2 = hm.get_point(c_bot, vclosest=[-601.99, 1483.13])
[c1] = hm.extract_subcontours(c_bot, [c[16], c[15]])
[cr, ct, cl] = hm.extract_subcontours(
    c1, [c[16], p1, p2, c[15]], "vertex")
[cb] = hm.extract_subcontours(g1, [c[16], c[15]], project_to="line")
ali = hm.connect_subcontours([cb, cr, ct, cl], fix=[0, 2], close="yes")

go = hm.add_unf_rect_grid(nx=3, ny=10)
g6 = hm.map_grid(
    go, ali,
    [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0.5]],
    [c[15], c[16], p1, p2, [-590, 1477]],
    project_to="vertex", algo="inverse_laplace")

# right input: g3
p0 = hm.get_point(g1cont, eclosest=c[18])
p1 = hm.get_point(c_bot, vclosest=[1785.79, 1489.68])
p2 = hm.get_point(c_bot, vclosest=[1776.04, 1489.65])
p3 = hm.get_point(g1cont, eclosest=c[17])
p4 = hm.get_point(g1cont, eclosest=[1748.54, 1481.69])
p5 = hm.get_point(g1cont, eclosest=[1742, 1482])
[c0, c1, c2, c4] = hm.extract_subcontours(
    c_bot, [p0, p1, p2, p3, p4], "line")
c3 = hm.create_contour([p0, p3])
c5 = hm.create_contour([p4, p5])
[c6] = hm.copy_geom(c4)
hm.scale_geom([c6], 150, 250)
ar1 = hm.connect_subcontours([c3, c0, c1, c2], fix=[0, 2], close="yes")
ar2 = hm.connect_subcontours([c3, c4, c5, c6], fix=[0, 2], close="yes")
go = hm.add_unf_rect_grid(nx=3, ny=10)
g31 = hm.map_grid(
    go, ar1,
    [[0, 0], [1, 0], [1, 1], [0, 1], [1, 0.5]],
    [p3, p0, p1, p2, [1774.22, 1483.7]],
    project_to="vertex", algo="inverse_laplace")
go = hm.add_unf_rect_grid(nx=3, ny=15)
g32 = hm.map_grid(
    go, ar2,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [p4, p5, p0, p3],
    project_to="vertex", algo="inverse_laplace")
g3 = hm.unite_grids(g31, [(g32, 0)])
# exclude superfluous zones in g1
[cout1] = hm.extract_subcontours(g3, [p2, p4])
cout2 = hm.create_contour([p2, p4])
cout = hm.unite_contours([cout1, cout2])
g1 = hm.exclude_contours(g1, cout, "inner")

# left top zone: g7
# 1. hole
u = hm.partition_contour(c_lhole, "const", step_bl)
bs = [0, step_bl / 5, 2 * step_bl / 5, 3 * step_bl / 5]
bgo = hm.BoundaryGridOptions(u, bs, "right", step_bl / 5, "no")
g71 = hm.build_boundary_grid(bgo)

# 2. sqr
c0 = hm.create_contour([c[8], c[13]])
c1 = hm.create_contour([c[13], c[19]])
c2 = hm.create_contour([c[21], c[22]])
c3 = hm.create_contour([c[20], c[23]])
c4 = hm.create_contour([c[19], c[32]])
[c5, c6, c81, c82, c83, c9, c10] = hm.extract_subcontours(
    c_ltop, [c[8], c[32], c[23], c[22], c[21],
             [117, 1924.7], [156, 1924.7], c[20]], "line")
c7 = hm.create_contour([c[19], c[20]])
c82 = hm.create_contour([c[22], c[21]])
c8 = hm.connect_subcontours([c81, c82, c83])

c0 = hm.partition_contour(c0, "const", step_bl, angle0=180, start=c[8])
nc0 = hm.info_contour(c0)['Nedges']
c1 = hm.partition_contour(c1, "const", step_bl, angle0=180)
nc1 = hm.info_contour(c1)['Nedges']
c3 = hm.partition_contour(c3, "const", step_bl, angle0=180, start=c[20],
                          nedges=nc0)
c4 = hm.partition_contour(c4, "const", step_bl, angle0=180, start=c[19],
                          nedges=nc0)
c5 = hm.partition_contour(c5, "const", 1, angle0=180, start=c[32],
                          nedges=nc1)
c6 = hm.partition_contour(c6, "const", 1, angle0=180, nedges=3)
c7 = hm.partition_contour(c7, "const", 1, angle0=180, nedges=3)
c8 = hm.partition_contour(c8, "const", step_bl, angle0=180, nedges=nc1)
c9 = hm.partition_contour(c9, "const", step_bl, angle0=180, nedges=nc0)
c10 = hm.partition_contour(c10, "const", step_bl, angle0=180, nedges=nc1)

g721 = hm.add_custom_rect_grid("inverse_laplace", c4, c1, c0, c5)
g722 = hm.add_custom_rect_grid("inverse_laplace", c3, c7, c4, c6)
g723 = hm.add_custom_rect_grid("inverse_laplace", c9, c10, c3, c8)

# boundary grid after sqr bottom side
bgo = hm.BoundaryGridOptions(c1, [0], "left", step_bl, "no")
bgo.uniform_partition(10 * step_bl, 10)
g73 = hm.build_boundary_grid(bgo)
mincont = hm.create_contour(
    [[163.373, 1962.943], [164.354, 1954.359],
     [177.761, 1945.612], [194.275, 1944.95],
     [195.010, 1950.517], [192.149, 1959.91],
     [212.823, 1989.198], [163.373, 1962.943]])
g73 = hm.exclude_contours(g73, mincont, "outer")
# 3. cavern
p1 = hm.get_point(g721, vclosest=[179, 1963])
[sa1] = hm.extract_subcontours(c_ltop, [c[19], c[13]])
sa1 = hm.partition_contour(sa1, "const", step_bl)
[sa2] = hm.extract_subcontours(g721, [c[19], p1])
[sa3] = hm.extract_subcontours(g721, [c[13], c[8]])
sa4 = hm.create_contour([c[8], p1])
sa1 = hm.unite_contours([sa1, sa2, sa3, sa4])

bs = [0, step_bl / 2, 2. * step_bl / 2, 3. * step_bl / 2, 4. * step_bl / 2]
bgo = hm.BoundaryGridOptions(sa1, bs, "left", step_bl, "no",
                             start_point=c[19], end_point=c[13])
g74 = hm.build_boundary_grid(bgo)

# 4. input
[cl, cb, cr, ct] = hm.extract_subcontours(c_ltop, [[105, 1968], [103, 1962],
                                                   c[21], c[22], [105, 1968]])
[cr] = hm.extract_subcontours(g723, [c[22], c[21]], "line")
nic = hm.connect_subcontours([cb, cl, ct, cr], [1, 3], "yes")
cl = hm.partition_contour(cl, "const", 1, nedges=3)
cb = hm.partition_contour(cb, "const", 1, nedges=3)
cr = hm.partition_contour(cr, "const", 1, nedges=3)
ct = hm.partition_contour(ct, "const", 1, nedges=3)
g75 = hm.add_custom_rect_grid("linear", cr, cb, cl, ct)

# 5. unite all
g7 = hm.unite_grids(g721, [(g722, 0), (g723, 0), (g73, 0)], buffer_fill='4')
g7 = hm.unite_grids(g7, [(g74, 1), (g75, 1)], buffer_fill='4')
g7 = hm.unite_grids(g7, [(g71, 1.2)], empty_holes=True, buffer_fill='4')

# right top zone: g8
# 1. hole
u = hm.partition_contour(c_rhole, "const", step_bl)
bs = [0, step_bl / 5, 2 * step_bl / 5, 3 * step_bl / 5]
bgo = hm.BoundaryGridOptions(u, bs, "right", step_bl / 5, "no")
g81 = hm.build_boundary_grid(bgo)

# 2. cavern
ccut = hm.create_contour([c[13], c[8], [137.5, 2011], c[32], c[19],
                          [158, 1935], [203, 1937], c[13]])
gorig = hm.exclude_contours(g7, ccut, "outer")
[c1] = hm.extract_subcontours(c_rtop, [c[34], c[7]])
c2 = hm.create_contour([c[7], c[14]])
[c3] = hm.extract_subcontours(c_rtop, [c[14], c[33]])
c4 = hm.create_contour([c[33], c[34]])
ctar = hm.unite_contours([c1, c2, c3, c4])
g82 = hm.map_grid(
    gorig, ctar,
    [c[13], c[8], c[32], [163, 1953]],
    [c[14], c[7], c[34], c[33]],
    project_to="vertex", is_reversed=True)

# 3. rest
[cl1] = hm.extract_subcontours(g82, [c[33], c[34]])
cl2 = hm.create_contour([c[33], c[37]])
cl2 = hm.partition_contour(cl2, "const", step_bl)
cl = hm.unite_contours([cl1, cl2])
cb = hm.create_contour([c[37], c[38]])
[cr, ct1, ct2, ct3] = hm.extract_subcontours(
    c_rtop, [c[38], c[39], c[36], c[35], c[34]])
ct2 = hm.create_contour([c[36], c[35]])
ct1 = hm.partition_contour(ct1, "const", step_bl)
ct2 = hm.partition_contour(ct2, "const", nedges=3)
ct3 = hm.partition_contour(ct3, "const", step_bl)
ct = hm.unite_contours([ct1, ct2, ct3])
cr = hm.partition_contour(cr, "const", step_bl,
                          nedges=hm.info_contour(cl)['Nedges'])
cb = hm.partition_contour(cb, "const", step_bl,
                          nedges=hm.info_contour(ct)['Nedges'])
g83 = hm.add_custom_rect_grid("inverse_laplace", cl, cb, cr, ct)

# 4. input
cl = hm.create_contour([c[36], c[35]])
[cout] = hm.extract_subcontours(c_rtop, [c[36], c[35]])
cout = hm.unite_contours([cl, cout])
go = hm.add_unf_rect_grid(nx=3, ny=3)
g84 = hm.map_grid(
    go, cout,
    [[0, 0], [1, 0], [1, 1], [0, 1]],
    [c[36], [1291.48, 1965.2], [1289.67, 1971.23], c[35]])

# 5. unite
g8 = hm.unite_grids(g82, [(g83, 0), (g84, 0)])
g8 = hm.unite_grids(g8, [(g81, 1)], empty_holes=True, buffer_fill='4',
                    zero_angle_approx=30)

# left-right top connect: g11
[cl, cb, cr, ct] = hm.extract_subcontours(
    c_top, [c[i] for i in [8, 13, 14, 7, 8]])
[cl] = hm.extract_subcontours(g7, [c[13], c[8]])
[cr] = hm.extract_subcontours(g8, [c[7], c[14]])
ct = hm.partition_contour(ct, "ref_weights", bstep, start=c[7])
cb = hm.partition_contour(cb, "ref_weights", bstep, start=c[14],
                          nedges=hm.info_contour(ct)['Nedges'])
g11 = hm.add_custom_rect_grid("inverse_laplace", cl, cb, cr, ct)

# top grid: g12
[bc, rc, tc, lc] = hm.extract_subcontours(
    c_mid, [c[i] for i in [5, 6, 14, 13, 5]])
[tc] = hm.extract_subcontours(g11, [c[13], c[14]])
lc = hm.partition_contour(lc, "ref_weights", bstep1, 180, start=c[13],
                          keep_pts=[c[9], c[10]])
rc = hm.partition_contour(rc, "ref_weights", bstep1, 180, start=c[14],
                          nedges=hm.info_contour(lc)['Nedges'],
                          keep_pts=[c[11], c[12]])
bc = hm.partition_contour(bc, "ref_weights", bstep, 180, start=c[5],
                          nedges=hm.info_contour(tc)['Nedges'])
g12 = hm.add_custom_rect_grid("inverse_laplace", lc, bc, rc, tc)

# left chairs: g9
[v1, h1, v2, ch1, v3, h2, v4, ch2] = hm.extract_subcontours(
    c_lchairs, [c[i] for i in [24, 28, 29, 25, 26, 30, 31, 27, 24]])
v1 = hm.partition_contour(v1, "ref_weights", bstep1, 180, start=c[28])
ne1 = hm.info_contour(v1)['Nedges']
v2 = hm.partition_contour(v2, "ref_weights", bstep1, 180, start=c[29],
                          nedges=ne1)
v3 = hm.partition_contour(v3, "ref_weights", bstep1, 180, start=c[30],
                          nedges=ne1)
v4 = hm.partition_contour(v4, "ref_weights", bstep1, 180, start=c[31],
                          nedges=ne1)
ch1 = hm.partition_contour(ch1, "const", step_main)
ch2 = hm.partition_contour(ch2, "const", step_main)
# main boundary layer
bcont = hm.connect_subcontours([v1, h1, v2, ch1, v3, h2, v4, ch2])
bgo1 = hm.BoundaryGridOptions(bcont, direction="right", bnd_stepping="no",
                              range_angles=[40, 125, 360, 360],
                              start_point=c[28], end_point=c[31])
bgo2 = hm.BoundaryGridOptions(bcont, direction="right", bnd_stepping="no",
                              range_angles=[40, 125, 360, 360],
                              start_point=c[30], end_point=c[29])
bgo1.incremental_partition(step_bl, 1.2, 4)
bgo2.incremental_partition(step_bl, 1.2, 4)
g911 = hm.build_boundary_grid([bgo1])
g912 = hm.build_boundary_grid([bgo2])

# fillers
c1 = hm.grid_bnd_to_contour(g911, False)
c2 = hm.grid_bnd_to_contour(g912, False)
[c1] = hm.extract_subcontours(c1, [c[28], c[31]])
[c2] = hm.extract_subcontours(c2, [c[30], c[29]])
[line1] = hm.extract_subcontours(c1, [c[29], c[24]])
[line2] = hm.extract_subcontours(c2, [c[28], c[25]])
[line3] = hm.extract_subcontours(c2, [c[31], c[26]])
[line4] = hm.extract_subcontours(c1, [c[30], c[27]])
line12b = hm.create_contour(
    [hm.get_point(line1, vclosest=c[28]),
     hm.get_point(line2, vclosest=c[29])])
line12t = hm.create_contour(
    [hm.get_point(line1, vclosest=c[24]),
     hm.get_point(line2, vclosest=c[25])])
line34b = hm.create_contour(
    [hm.get_point(line3, vclosest=c[30]),
     hm.get_point(line4, vclosest=c[31])])
line34t = hm.create_contour(
    [hm.get_point(line3, vclosest=c[26]),
     hm.get_point(line4, vclosest=c[27])])
line12t = hm.partition_contour(line12t, "const", 1.3 * step_main)
ne = hm.info_contour(line12t)['Nedges']
line12b = hm.partition_contour(line12b, "const", nedges=ne)
line34t = hm.partition_contour(line34t, "const", nedges=ne)
line34b = hm.partition_contour(line34b, "const", nedges=ne)
g921 = hm.add_custom_rect_grid(
    "inverse_laplace", line1, line12b, line2, line12t)
g922 = hm.add_custom_rect_grid(
    "inverse_laplace", line3, line34b, line4, line34t)

g9 = hm.unite_grids(g911, [(g912, 0), (g921, 0), (g922, 0)])

# right chairs: g10
[v1, h1, v2, ch1, b1, b2, b3, b4, b5, ch2] = hm.extract_subcontours(
    c_rchairs, [c[i] for i in [40, 42, 43, 41, 44, 45, 46, 48, 47, 49, 40]])
p1 = hm.get_point(g1cont, eclosest=c[44])
p2 = hm.get_point(g1cont, eclosest=c[46])
p3 = hm.get_point(g1cont, eclosest=c[48])
p4 = hm.get_point(g1cont, eclosest=c[49])
[cg1, cg2, cg3] = hm.extract_subcontours(g1cont, [p1, p2, p3, p4], "line")
c1 = hm.create_contour([p1, c[44]])
c2 = hm.create_contour([p2, c[46]])
c3 = hm.create_contour([p3, c[48]])
c4 = hm.create_contour([p4, c[49]])

v1 = hm.partition_contour(v1, "ref_weights", bstep1, 180, start=c[42])
ne1 = hm.info_contour(v1)['Nedges']
v2 = hm.partition_contour(v2, "ref_weights", bstep1, 180, start=c[43],
                          nedges=ne1)
ch1 = hm.partition_contour(ch1, "const", step_main)
ch2 = hm.partition_contour(ch2, "const", step_main)
b1 = hm.partition_contour(b1, "const", step_main)
b3 = hm.partition_contour(b3, "const", step_main)
b5 = hm.partition_contour(b5, "const", step_main)
cg1 = hm.partition_contour(cg1, "const", nedges=hm.info_contour(b1)['Nedges'])
cg2 = hm.partition_contour(cg2, "const", nedges=hm.info_contour(b3)['Nedges'])
cg3 = hm.partition_contour(cg3, "const", nedges=hm.info_contour(b5)['Nedges'])

# 1. right zone
c1 = hm.partition_contour(c1, "ref_points",
                          [4 * step_bl, p1, 5 * step_bl, c[44]])
nc1 = hm.info_contour(c1)['Nedges']
c2 = hm.partition_contour(c2, "const", 4 * step_bl)
nc2 = hm.info_contour(c2)['Nedges']
b2 = hm.partition_contour(b2, "const", nedges=nc1 - nc2)
cb = hm.unite_contours([b2, c2])
g101 = hm.add_custom_rect_grid("inverse_laplace", b1, cb, cg1, c1)
c3 = hm.partition_contour(c3, "const", nedges=hm.info_contour(c2)['Nedges'])
g102 = hm.add_custom_rect_grid("inverse_laplace", c3, cg2, c2, b3)
c4 = hm.partition_contour(c4, "const", nedges=hm.info_contour(c3)['Nedges'])
pm = hm.get_point(c_rchairs, eclosest=[2147.5, 396.3])
cm = hm.create_contour([c[48], pm, c[49]])
cm = hm.partition_contour(cm, "const", nedges=hm.info_contour(cg3)['Nedges'])
g103 = hm.add_custom_rect_grid("inverse_laplace", cm, c4, cg3, c3)
g104 = hm.add_triangle_grid(c[48], c[47], pm, 3)

# 2. main boundary layer
[cupper, _, clower] = hm.extract_subcontours(
    g1cont, [c[2], p1, p4, c[1]], "line")
bcont = hm.connect_subcontours([clower, c4, ch2, v1, h1, v2, ch1, c1, cupper])
bgo1 = hm.BoundaryGridOptions(bcont, direction="right", bnd_stepping="no",
                              range_angles=[40, 95, 360, 360],
                              start_point=p1, end_point=c[43])
bgo2 = hm.BoundaryGridOptions(bcont, direction="right", bnd_stepping="no",
                              range_angles=[40, 95, 360, 360],
                              start_point=c[42], end_point=p4)
bgo1.incremental_partition(step_bl, 1.2, 4)
bgo2.incremental_partition(step_bl, 1.2, 4)
g105 = hm.build_boundary_grid([bgo1])
g106 = hm.build_boundary_grid([bgo2])

# 3. fillers
c1 = hm.grid_bnd_to_contour(g105, False)
c2 = hm.grid_bnd_to_contour(g106, False)
[c1] = hm.extract_subcontours(c1, [c[44], c[42]])
[c2] = hm.extract_subcontours(c2, [c[43], c[49]])
[line1] = hm.extract_subcontours(c2, [c[43], c[40]])
[line2] = hm.extract_subcontours(c1, [c[41], c[42]])
line12b = hm.create_contour(
    [hm.get_point(line1, vclosest=c[42]),
     hm.get_point(line2, vclosest=c[43])])
line12t = hm.create_contour(
    [hm.get_point(line1, vclosest=c[40]),
     hm.get_point(line2, vclosest=c[41])])
line12t = hm.partition_contour(line12t, "const", 1.3 * step_main)
ne = hm.info_contour(line12t)['Nedges']
line12b = hm.partition_contour(line12b, "const", nedges=ne)
g107 = hm.add_custom_rect_grid(
    "inverse_laplace", line1, line12b, line2, line12t)

g10 = hm.unite_grids(
    g101, [(g, 0) for g in [g102, g103, g104, g105, g106, g107]])

excont = hm.create_contour([c[44], p1, [2423, 734], p4, c[49], c[44]])
g1 = hm.exclude_contours(g1, excont, "inner")

# assemble all
g = g1
print "UNITE 1. chairs"
g = hm.unite_grids(g, [(g9, 10)], zero_angle_approx=10,
                   empty_holes=True, buffer_fill='4')
g = hm.unite_grids(g, [(g10, 10)], zero_angle_approx=10,
                   empty_holes=True, buffer_fill='4')
print "UNITE 2. input"
g = hm.unite_grids(g, [(g6, 10)], zero_angle_approx=10, buffer_fill='4')
g = hm.unite_grids(g, [(g3, 10)], zero_angle_approx=10, buffer_fill='4')
print "UNITE 3. middle and top"
g = hm.unite_grids(g, [(g12, 20)], zero_angle_approx=10, buffer_fill='4')
g = hm.unite_grids(g, [(g11, 0)])
print "UNITE 4. appendices"
g = hm.unite_grids(g, [(g4, 5)], zero_angle_approx=10, buffer_fill='4')
g = hm.unite_grids(g, [(g5, 5)], zero_angle_approx=10, buffer_fill='4')
print "UNITE 5. top zones"
g = hm.unite_grids(g, [(g7, 0)])
g = hm.unite_grids(g, [(g8, 0)])

hm.export_grid_vtk([g], "a.vtk")
