from hybmeshpack import hmscript as hm
hm.check_compatibility('0.4.6', 2)

# boundary features
brib = hm.add_boundary_type(1, "rib")
bbottom = hm.add_boundary_type(2, "bottom")
binput = hm.add_boundary_type(3, "input")
boutput = hm.add_boundary_type(4, "output")
btop = hm.add_boundary_type(5, "top")

# ===================== Options
# area sizes
Lx1 = 16.66
Lx2 = 25.
Ly = 3.33
Lz = 8.33
hcyl = 0.5

# auxiliary sizes
Lx1_rib_blay = 0.15
Lx2_rib_blay = 0.7
Ly_rib_blay = 0.3
Lx_wake_reg = 15.
Ly_wake_reg = 1.0

# step sizes
step0 = 0.03  # lowest step in the walls boundary layers
step1 = 0.03  # step in the boudary layer straight behind the rib
step2 = 0.02 # outer rib boundary region step
step3 = 0.15 # maximum wake region step, maximum wall boundary step
step4 = 0.25  # largest step

# ======================= Meshing procedure
# 1. channel region
region_channel = hm.add_rect_contour([-Lx1, 0], [Lx2, Ly],
                                     [bbottom, boutput, btop, binput])

# 2. rib boundary layer region
cyl = hm.add_circ_contour([0, 0], hcyl / 2, 128, brib)
[bcyl] = hm.copy_geom([cyl])
xpc = 100. * (hcyl + Lx1_rib_blay + Lx2_rib_blay) / (hcyl)
ypc = 100. * (hcyl + 2 * Ly_rib_blay) / (hcyl)
hm.scale_geom(bcyl, xpc, ypc)
xmv = (Lx2_rib_blay - Lx1_rib_blay) / 2.
hm.move_geom(bcyl, xmv, 0.)
region_blay = hm.unite_contours([cyl, bcyl])
region_blay = hm.clip_domain(region_blay, region_channel, "intersection")

# 3. wake region
x0 = 1.1 * hcyl / 2
x1 = Lx_wake_reg
y0 = 0
y10 = 0.9 * (Ly_rib_blay + hcyl / 2)
y11 = Ly_wake_reg
region_wake = hm.create_contour([[x0, y0], [x1, y0], [x1, y11], [x0, y10],
                                 [x0, y0]])
region_wake = hm.clip_domain(region_wake, bcyl, "difference")
pcross = hm.get_point(region_wake, vclosest=[-Lx1, Ly])

# 4. meshing bcyl region
p1 = hm.get_point(region_blay, vclosest=[-Lx1, 0])
p2 = hm.get_point(region_blay, vclosest=[-hcyl / 2, 0])
p3 = hm.get_point(region_blay, vclosest=[hcyl / 2, 0])
p4 = hm.get_point(region_blay, vclosest=[Lx2, 0])
[b1, b2, b3, b4] = hm.extract_subcontours(region_blay, [p1, p2, p3, p4, p1])
b3 = hm.partition_contour(
    b3, "ref_weights", [step1, 0, step2, 1], start=p3)
b4 = hm.partition_contour(
    b4, "ref_points", [step2, p4, step2, pcross, step2, p1],
    keep_pts=[pcross])
mesh_blay = hm.add_custom_rect_grid("orthogonal", b3, b4, b1, b2)

# 5. meshing wake region
# 5.1 separation into transitional traingle and main region
p1 = [hcyl / 2 + Lx2_rib_blay, 0]
p2 = [p1[0] + hcyl / 2, Ly_wake_reg]
ctmp = hm.create_contour([p1, p2])
region_wake = hm.unite_contours([region_wake, ctmp])
rr = hm.decompose_contour(region_wake)
region_wake_trans = hm.pick_contour(pcross, rr)
region_wake_main = hm.pick_contour([Lx2, 0], rr)
p2 = hm.get_point(region_wake_trans, vclosest=[Lx2, Ly])
[b1, b2, b3] = hm.extract_subcontours(region_wake_trans, [p1, p2, pcross, p1])
[b3] = hm.extract_subcontours(mesh_blay, [p1, pcross])
region_wake_trans = hm.connect_subcontours([b1, b2, b3])
# 5.2 meshing transitional triangle
[b3] = hm.extract_subcontours(region_wake_trans, [pcross, p1])
nedges = hm.info_contour(b3)['Nedges']
b3pts = [hm.get_point(b3, i) for i in range(nedges + 1)]
b3pts.sort(key=lambda a: a[1])
b3pts_source = [[0., float(i) / nedges] for i in range(nedges + 1)]
b3pts.append(p2)
b3pts_source.append([1., 0.])
tri = hm.add_triangle_grid([0, 0], [1, 0], [0, 1], nedges)
mesh_wake_trans = hm.map_grid(
    tri, region_wake_trans, b3pts_source, b3pts, snap="shift_vertices")
# 5.3 meshing main wake region
[b1, b2, b3, b4] = hm.extract_subcontours(
    region_wake_main,
    [p1, [Lx_wake_reg, 0.], [Lx_wake_reg, Ly_wake_reg], p2, p1])
[b4] = hm.extract_subcontours(mesh_wake_trans, [p1, p2])
b1 = hm.partition_contour(b1, "ref_weights", [step2, 0, step3, 1], start=p1)
mesh_wake_main = hm.add_custom_rect_grid("orthogonal", b4, b1, b2, b3)
# 5.4 connect all
mesh_around = hm.unite_grids(mesh_wake_main, [(mesh_wake_trans, 0),
                                              (mesh_blay, 0)])
hm.reflect_geom([mesh_wake_trans, mesh_blay], [0, 0], [1, 0])
mesh_around = hm.unite_grids(mesh_around, [(mesh_wake_trans, 0),
                                           (mesh_blay, 0)])
# 5.5 bottom wake part
cmesh_around = hm.grid_bnd_to_contour(mesh_around)
p1 = hm.get_point(cmesh_around, vclosest=[hcyl / 2 + Lx2_rib_blay, -Ly])
p2 = [Lx_wake_reg, p1[1]]
p3 = [p2[0], 0]
p4 = [hcyl / 2 + Lx2_rib_blay, 0]
[bl, bt] = hm.extract_subcontours(mesh_around, [p1, p4, p3])
bb = hm.create_contour([p1, p2])
br = hm.create_contour([p2, p3])
mesh_bot_wake = hm.add_custom_rect_grid("orthogonal", bt, bl, bb, br)
mesh_around = hm.unite_grids(mesh_around, [(mesh_bot_wake, 0)])

if not hm.skewness(mesh_around, 0.8)['ok']:
    raise Exception

# # 5.6 boundary layer under the wake
# p4 = hm.get_point(cmesh_around, vclosest=[-hcyl / 2, -Ly_rib_blay - hcyl / 2])
# p3 = p2
# p1 = [p4[0], -1.5 * hcyl]
# p2 = [Lx_wake_reg, p1[1]]
# [bt] = hm.extract_subcontours(mesh_around, [p4, p3])
# bl = hm.create_contour([p1, p4])
# bb = hm.create_contour([p1, p2])
# br = hm.create_contour([p2, p3])
# bl = hm.partition_contour(bl, "ref_weights", [step0, 0, step2, 1], start=p1)
# mesh_under = hm.add_custom_rect_grid("orthogonal", bt, bl, bb, br)
# mesh_around = hm.unite_grids(mesh_around, [(mesh_under, 0)])

# # 5.7 move up
# hm.move_geom(mesh_around, 0, 1.5 * hcyl)

# # 6. Wall boundary layer
# bh = 1.5 * hcyl
# xmid = 0
# bbot1 = hm.create_contour([[-Lx1, 0], [xmid, 0]])
# bbot2 = hm.create_contour([[xmid, 0], [Lx2, 0]])
# bleft = hm.create_contour([[-Lx1, 0], [-Lx1, step1]])
# bmid = hm.create_contour([[xmid, 0], [xmid, bh]])
# bright = hm.create_contour([[Lx2, 0], [Lx2, bh]])

# part_1 = [step0, 0, step2, hcyl / Lx1, step4, Lx1 / Lx1]
# part_2 = [step0, 0, step3, 1]

# bmid = hm.partition_contour(bmid, "ref_weights", part_2, start=[0, 0])
# ny = hm.info_contour(bmid)['Nedges']
# bleft = hm.partition_contour(bleft, "ref_weights", part_2, start=[-Lx1, 0],
#                              nedges=ny)
# bright = hm.partition_contour(bright, "ref_weights", part_2, start=[Lx2, 0],
#                               nedges=ny)
# bbot1 = hm.partition_contour(bbot1, "ref_weights", part_1, start=[-Lx1, 0])
# bbot2 = hm.partition_contour(bbot2, "const", step4)

# wallbb1 = hm.add_custom_rect_grid("linear", bbot1, bmid, None, bleft)
# wallbb2 = hm.add_custom_rect_grid("linear", bbot2, bright, None, bmid)
# wallbb3 = hm.unite_grids(wallbb1, [(wallbb2, 0)])
# [wallbb4] = hm.copy_geom([wallbb3])
# hm.reflect_geom([wallbb4], [0, Ly / 2], [1, Ly / 2])
# mesh_wall_blay = hm.unite_grids(wallbb3, [(wallbb4, 0)])

# # 7. Channel
# region_outer = hm.clip_domain(region_channel, mesh_wall_blay, "difference",
#                               simplify=False)
# [b1, b2, b3, b4] = hm.extract_subcontours(
#     region_outer, [[-Lx1, 0], [Lx2, 0], [Lx2, Ly], [-Lx1, Ly], [-Lx1, 0]])
# b2 = hm.partition_contour(b2, "const", step4)
# b4 = hm.partition_contour(b4, "const", nedges=hm.info_contour(b2)['Nedges'])
# mesh_outer = hm.add_custom_rect_grid("direct_laplace", b1, b2, b3, b4)

# # 8. Superposition
# g1 = hm.unite_grids(mesh_outer, [(mesh_wall_blay, 0)])
# final = hm.unite_grids(g1, [(mesh_around, 0.95 * step4)], buffer_fill='4')
# outreg = hm.decompose_contour(mesh_around)
# outreg = hm.pick_contour([0, 1.5 * hcyl], outreg)
# final = hm.exclude_contours(final, outreg, "inner")
# final = hm.exclude_contours(final, region_channel, "outer")
# hm.heal_grid(final)

# # 9. Export
# hm.export_grid_vtk([final], "fullcyl.vtk")
# hm.export_contour_vtk([final], "fullcyl_cont.vtk")

# # 10. Checks
# print "Skewness", hm.skewness(final)
# print "Number of cells:", hm.info_grid(final)['Ncells']
