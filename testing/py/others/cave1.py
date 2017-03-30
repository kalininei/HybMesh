from hybmeshpack import hmscript as hm
from math import *  # NOQA
hm.check_compatibility('0.4.6', 2)
global x2, y2

#area = hm.import_contour_hmc('cave1.hmc')
#hm.export_contour_vtk(area, "init.vtk")

#p1 = hm.get_point(area, vclosest=[227, 0])
#p2 = hm.get_point(area, vclosest=[278, 0])
#p3 = hm.get_point(area, vclosest=[252, -50])
#a1 = [0, 0]
#a2 = [500, 0]
#a3 = [500, 50]
#a4 = [0, 50]

bcirc = hm.add_boundary_type(1, "cave")
binput = hm.add_boundary_type(2, "input")
boutput = hm.add_boundary_type(3, "outflow")
bbot = hm.add_boundary_type(4, "bottom")
btop = hm.add_boundary_type(5, "top")

# option 1
# circ = hm.add_circ_contour2(p1, p2, p3, 128, bcirc)
# [a, b] = hm.extract_subcontours(circ, [p1, p2, p1], "line")
# circ = hm.connect_subcontours([a, b])
# p1 = hm.get_point(circ, vclosest=p1)
# p2 = hm.get_point(circ, vclosest=p2)
# hm.export_contour_vtk(circ, "c1.vtk")

# grbase = hm.add_circ_rect_grid([0, 0], 1., 0.1, 1.0, 0.1, "orthogonal_circ")
# x = math.sqrt(1. - 0.5 * 0.5)
# grmapped = hm.map_grid(grbase, circ, [[-x, 0.5], [x, 0.5]], [p1, p2],
#                        project_to="vertex", btypes="from_contour")
# hi = hm.info_contour(grmapped)
# h = hi['length'] / hi['Nedges']

# b1 = hm.create_contour([a1, p1])
# [b2] = hm.extract_subcontours(grmapped, [p2, p1])
# b3 = hm.create_contour([p2, a2])
# b1 = hm.partition_contour(b1, "const", h)
# b3 = hm.partition_contour(b3, "const", h)
# l = hm.create_contour([a1, a4], binput)
# r = hm.create_contour([a2, a3], boutput)
# t = hm.create_contour([a4, a3], btop)
# lpart = [h / 10., 0., h, 20., h, -20., h / 10., -1e-10]
# l = hm.partition_contour(l, "ref_lengths", lpart, start=a1)
# b = hm.connect_subcontours([b1, b2, b3])
# hm.set_boundary_type(b, bbot)
# mesh_area = hm.add_custom_rect_grid("orthogonal", l, b, r, t)
# case1 = hm.unite_grids(grmapped, [(mesh_area, 0)])

# hm.export_grid_vtk([case1], "cave1_1.vtk")
# hm.export_contour_vtk([case1], "cave1_1_cont.vtk")
# hm.export_grid_msh(case1, "cave1_1.msh")

# option 2
# ar = hm.create_contour([a1, a2, a3, a4, a1])
# hm.set_boundary_type(ar, [bbot, boutput, btop, binput])
# ar = hm.clip_domain(ar, circ, "union")
# [b1] = hm.extract_subcontours(ar, [a1, a2])
# [b2] = hm.extract_subcontours(ar, [a3, a4])
# b1 = hm.partition_contour(b1, "const", h)
# b2 = hm.partition_contour(b2, "const", h)
# bopt = hm.BoundaryGridOptions(b1, bnd_stepping="no")
# bopt.incremental_partition(h / 10, 1.11, 15)
# bopt.partition.append(bopt.partition[-1] + 0.7 * h)
# g1 = hm.build_boundary_grid(bopt)
# bopt.contour_id = b2
# g2 = hm.build_boundary_grid(bopt)
# g3 = hm.add_unf_rect_grid([a1[0], -a4[1]], a3, custom_x=h, custom_y=h)

# [b1] = hm.extract_subcontours(g1, [a1, a2])
# [b2] = hm.extract_subcontours(g2, [a3, a4])
# l = hm.create_contour([a1, a4], binput)
# r = hm.create_contour([a2, a3], boutput)
# ar2 = hm.connect_subcontours([b1, r, b2, l])

# g3 = hm.exclude_contours(g3, ar2, "outer")
# g1 = hm.exclude_contours(g1, ar2, "outer")
# g2 = hm.exclude_contours(g2, ar2, "outer")

# case2 = hm.unite_grids(g3, [(g1, h / 2.), (g2, h / 2.)], buffer_fill='4')

# hm.export_grid_vtk([case2], "cave1_2.vtk")
# hm.export_contour_vtk([case2], "cave1_2_cont.vtk")
# hm.export_grid_msh(case2, "cave1_2.msh")

# option 3
# geometry
p1 = [226.1, 0]
p2 = [278.262, 0]
rad = 30
pc = [252.182, -15.1]
p1[0] = -sqrt(rad*rad - (p1[1]-pc[1])*(p1[1]-pc[1])) + pc[0]
p2[0] = sqrt(rad*rad - (p2[1]-pc[1])*(p2[1]-pc[1])) + pc[0]
alpha = atan2(p2[1]-pc[1], p2[0]-pc[0])
x2 = 500
y2 = 50

# parameters
# cavern parameters
#hmin = 0.2
#hcave = 2.0
#shiftlay = 2.0
hmin = 0.4
hhormin = 3*hmin
hblayer = 1.5
hcave = 4.0
shiftlay = 3.0
layheight = 7
blay = hm.partition_segment(0, layheight, hmin, hblayer)
buf = 2
internal_step = 1.5*hblayer

# main area parameters
hmain = 5.

arc = hm.partition_segment(0, alpha, hblayer/rad, hhormin/rad) +\
      hm.partition_segment(alpha, pi-alpha, hhormin/rad, hhormin/rad,
                           [pi/2., hcave/rad])[1:]+\
      hm.partition_segment(pi-alpha, pi, hhormin/rad, hblayer/rad)[1:]+\
      hm.partition_segment(pi, 2.*pi, hblayer/rad, hblayer/rad,
                           [1.5*pi, hcave/rad])[1:]

rads = [0, rad-layheight]
for b in list(reversed(blay))[1:]:
    rads.append(rad-b)

g = hm.add_unf_circ_grid(pc, custom_arcs=arc, custom_rads=rads,
                         is_trian=False)

gfiller = hm.add_circ_rect_grid(pc, rad-1.3*layheight, internal_step)
g = hm.unite_grids(g, [(gfiller, buf)])

chord = hm.add_circ_contour2(p2, [pc[0], p1[1] + shiftlay], p1, 512)
mapping_area = hm.clip_domain(chord, g, "intersection")
hm.export_contour_vtk(mapping_area, "c1.vtk")

g2 = hm.map_grid(g, mapping_area, [p1, p2], [p1, p2], algo="direct_laplace")

pbef, paft = [p1[0]-layheight, 0], [p2[0]+layheight, 0]
c1 = hm.create_contour([pbef, p1])
c1 = hm.partition_contour(c1, "ref_weights", [hblayer, 1, hhormin, 0], start=p1)
c2 = hm.create_contour([p2, paft])
c2 = hm.partition_contour(c2, "ref_weights", [hhormin, 0, hblayer, 1], start=p2)
[c3] = hm.extract_subcontours(g2, [p2, p1])
c4 = hm.create_contour([pbef, [pbef[0], 100]])
c5 = hm.create_contour([[paft[0], 100], paft])
c3 = hm.connect_subcontours([c4, c1, c3, c2, c5], [2], shiftnext=False)
[c3] = hm.decompose_contour(c3)
bopt = hm.BoundaryGridOptions(c3, blay, bnd_stepping="no",
    start_point=pbef, end_point=paft)
g3 = hm.build_boundary_grid([bopt])
hm.set_boundary_type(g2, 1)
hm.set_boundary_type(g3, 4)
g2 = hm.unite_grids(g2, [(g3, 0)])

# main area
bot = hm.create_contour([[0, 0], pbef, paft, [x2, 0]])
bot = hm.partition_contour(bot, algo="const", step=hmain, angle0=-1)
left = map(
    lambda x: [0, x],
    blay[:-1] + hm.partition_segment(layheight, y2,
                                     hblayer, hmin,
                                     [y2/2., hmain, y2-layheight, hblayer]))
left = hm.create_contour(left)
outh = hm.get_point(left, vclosest=[0, paft[1]+layheight+shiftlay+2*hmain])[1]
[right] = hm.copy_geom(left)
[top] = hm.copy_geom(bot)
hm.set_boundary_type(left, 2)
hm.set_boundary_type(right, 3)
hm.set_boundary_type(bot, 4)
hm.set_boundary_type(top, 5)
gmain = hm.add_custom_rect_grid("linear", left, bot, right, top)

outcont = hm.add_rect_contour(pbef, [paft[0], outh])
gmain = hm.exclude_contours(gmain, outcont, "inner")

#final unite
gout = hm.unite_grids(gmain, [(g2, 0)])
hm.export_contour_vtk(gout, "c2.vtk")
cc = hm.simplify_contour(gout, simplify=False, separate=True)
filler = hm.pick_contour([pc[0], outh], cc)
gg = hm.triangulate_domain(filler)
gout = hm.unite_grids(gout, [(gg, 0)])

#zcoords = hm.partition_segment(0, 50, 3, 3)
#gm3 = hm.extrude_grid(gout, zcoords)
#hm.export3d_grid_vtk(gm3, "gm3.vtk")

hm.export_grid_vtk(gout, "cave1_3.vtk")
hm.export_contour_vtk(gout, "cave1_3_cont.vtk")
hm.export_grid_msh(gout, "cave1_3.msh")
