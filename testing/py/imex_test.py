# import os.path
from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript import _dbg as hmdbg
hm.check_compatibility("0.4.6")

print "export to natives"
# export/import grids 2d to native
g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 3, 4)
g2 = hm.add_unf_rect_grid([0, 0], [1, 1], 2, 3)
g3 = hm.add_unf_rect_grid([0, 0], [1, 1], 5, 6)
hm.export_grid_hmg(g1, "g1.hmg")
hm.export_grid_hmg(g1, "g2.hmg", "bin")
hm.export_grid_hmg(g1, "g3.hmg", "fbin")
hm.export_grid_hmg([g1, g2, g3], "g4.hmg")
g4 = hm.import_grid_hmg("g2.hmg")
hmdbg.check(hm.info_grid(g4)['Ncells'] == 12)
g6 = hm.import_grid_hmg("g4.hmg", allgrids=True)
hmdbg.check(len(g6) == 3 and hm.info_grid(g6[1])['Ncells'] == 6)
for i, z in enumerate(hm.registered_grids()):
    hmdbg.check(int(z[-1]) == i + 1)
hm.remove_all()
g7 = hm.import_grid_hmg("g4.hmg", "Grid3")
hmdbg.check(g7 == "Grid3" and len(hm.registered_grids()) == 1)
hmdbg.check(hm.info_grid(g7)['Ncells'] == 30)

# export/import contours 2d to native
c1 = hm.add_rect_contour([0, 0], [1, 1], 1)
c2 = hm.add_rect_contour([0, 0], [2, 2], 2)
hm.export_contour_hmc(c1, "c1.hmc")
hm.export_contour_hmc(c2, "c2.hmc", "bin")
hm.export_contour_hmc([c1, c2], "c3.hmc", "bin")
c3 = hm.import_contour_hmc("c1.hmc", "Contour1")
[c4, c5] = hm.import_contour_hmc("c3.hmc", allconts=True)
c6 = hm.import_contour_hmc("c3.hmc")
hmdbg.check(len(hm.registered_contours()) == 6)
hmdbg.check(1 in hm.info_contour(c3)['btypes'])
hmdbg.check(1 in hm.info_contour(c4)['btypes'])
hmdbg.check(2 in hm.info_contour(c5)['btypes'])
hmdbg.check(1 in hm.info_contour(c6)['btypes'])

# export/import grids 3d to native
gg1 = hm.extrude_grid(g7, [0, 1])
gg2 = hm.extrude_grid(g7, [0, 1, 2])
hm.export3d_grid_hmg(gg1, "gg1.hmg", "ascii")
hm.export3d_grid_hmg([gg1, gg2], "gg2.hmg", "bin")
[gg3] = hm.import3d_grid_hmg("gg1.hmg", allgrids=True)
[gg4, gg5] = hm.import3d_grid_hmg("gg2.hmg", allgrids=True)
gg6 = hm.import3d_grid_hmg("gg2.hmg", "Grid3D_2")
hmdbg.check(len(hm.registered_grids3d()) == 6)
hmdbg.check(hm.info_grid3d(gg3)['Ncells'] == 30)
hmdbg.check(hm.info_grid3d(gg4)['Ncells'] == 30)
hmdbg.check(hm.info_grid3d(gg5)['Ncells'] == 60)
hmdbg.check(hm.info_grid3d(gg6)['Ncells'] == 60)

# export/import surfaces 3d to native
[ss1] = hm.grid3_bnd_to_surface(gg1, separate=True)
ss2 = hm.grid3_bnd_to_surface(gg2)
hm.export3d_surface_hmc(ss1, "ss1.hmc", "ascii")
hm.export3d_surface_hmc([ss2, ss1], "ss2.hmc", "bin")
hm.export3d_surface_hmc([ss1, gg1, ss2], "ss3.hmc", "fbin")
[ss3, ss4] = hm.import3d_surface_hmc("ss2.hmc", allsurfs=True)
ss5 = hm.import3d_surface_hmc("ss2.hmc")
ss6 = hm.import3d_surface_hmc("ss3.hmc", "SurfaceOfGrid3D_1")
hmdbg.check(len(hm.registered_surfaces()) == 6)
hmdbg.check(hm.info_surface(ss1)['Nfaces'] == hm.info_surface(gg1)['Nfaces'])
hmdbg.check(hm.info_surface(ss2)['Nnodes'] == hm.info_surface(gg2)['Nnodes'])
hmdbg.check(hm.info_surface(ss3)['Nfaces'] == hm.info_surface(gg2)['Nfaces'])
hmdbg.check(hm.info_surface(ss4)['Nedges'] == hm.info_surface(gg1)['Nedges'])
hmdbg.check(hm.info_surface(ss5)['Nfaces'] == hm.info_surface(gg2)['Nfaces'])
hmdbg.check(hm.info_surface(ss6)['Nedges'] == hm.info_surface(gg1)['Nedges'])

# export/import all to native
zc2 = hm.registered_contours()
zg2 = hm.registered_grids()
zg3 = hm.registered_grids3d()
zs3 = hm.registered_surfaces()
hm.export_all_hmd("all.hmd", fmt="bin")
hm.remove_all()
hmdbg.check(not hm.registered_contours())
hmdbg.check(not hm.registered_grids())
hmdbg.check(not hm.registered_grids3d())
hmdbg.check(not hm.registered_surfaces())
a = hm.import_all_hmd("all.hmd")
hmdbg.check(len(a[0] + a[1] + a[3] + a[3]) == len(zc2 + zg2 + zg3 + zs3))
hmdbg.check(all(a == b for a, b in zip(zc2, hm.registered_contours())))
hmdbg.check(all(a == b for a, b in zip(zg2, hm.registered_grids())))
hmdbg.check(all(a == b for a, b in zip(zs3, hm.registered_surfaces())))
hmdbg.check(all(a == b for a, b in zip(zg3, hm.registered_grids3d())))

# save/load project
hm.add_boundary_type(15, "boundary15")
numcom = hm.flow.com_count()
hm.save_project("p.hmp", "bin")
hm.load_project("../external_files/empty.hmp")
hmdbg.check(not hm.registered_contours())
hmdbg.check(not hm.registered_grids())
hmdbg.check(not hm.registered_grids3d())
hmdbg.check(not hm.registered_surfaces())
hmdbg.check(hm.registered_btypes() == [(0, 'None')])

hm.load_project("p.hmp")
hmdbg.check(numcom == hm.flow.com_count())
hmdbg.check(all(a == b for a, b in zip(zc2, hm.registered_contours())))
hmdbg.check(all(a == b for a, b in zip(zg2, hm.registered_grids())))
hmdbg.check(all(a == b for a, b in zip(zs3, hm.registered_surfaces())))
hmdbg.check(all(a == b for a, b in zip(zg3, hm.registered_grids3d())))
hmdbg.check(hm.registered_btypes() == [(0, 'None'), (15, 'boundary15')])
hm.load_project("../external_files/empty.hmp")

print "export 2d to fluent"
bleft = hm.add_boundary_type(1, "bleft")
bright = hm.add_boundary_type(2, "bright")
bbot = hm.add_boundary_type(3, "bbot")
btop = hm.add_boundary_type(4, "btop")
g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10)
g2 = hm.add_triangle_grid([0.3, -0.4], [0.5, 0.2], [0.7, -0.4], 10)
[g3] = hm.copy_geom(g2)
hm.reflect_geom(g3, [0, 0.5], [1, 0.5])
g4 = hm.unite_grids(g1, [(g2, 0.05), (g3, 0.05)])
hm.set_boundary_type(g4, bfun=lambda x0, y0, x1, y1, b:
                     {True: 0,
                      x0 + x1 == 0: 1,
                      x0 + x1 == 2: 2,
                      y0 + y1 == -0.8: 3,
                      y0 + y1 == 2.8: 4}[True])
hm.export_grid_vtk(g4, "g1.vtk")
hm.export_contour_vtk(g4, "c1.vtk")
hmdbg.check_ascii_file(15697319238564148717, "g1.vtk", "dev")
hmdbg.check_ascii_file(16408920837426241157, "c1.vtk", "dev")
hm.export_grid_msh(g4, "g1.msh")
hmdbg.check_ascii_file(17685805227099775273, "g1.msh", "dev")

print "export 2d to fluent with periodic conditions"
hm.export_grid_msh(g4, "g1.msh", [bbot, btop, True])
hmdbg.check_ascii_file(3165188744058895474, "g1.msh", "dev")

print "export 2d to fluent with double periodic conditions"
hm.export_grid_msh(g4, "g1.msh", [bbot, btop, True, bleft, bright, True])
hmdbg.check_ascii_file(2759680993089544531, "g1.msh", "dev")

print "controlled fail on illegal periodic data"
try:
    hm.export_grid_msh(g4, "g1.msh", [bbot, 0, True])
    hmdbg.check(False)
except hm.ExportError:
    hmdbg.check(True)

print "export to gmsh file"
hm.export_grid_gmsh(g4, "g1.msh")
hmdbg.check_ascii_file(15544215273974397325, "g1.msh", "dev")

print "import from gmsh file"
hm.remove_all()
p1, p2, p3, p4, p5 = [0, 0], [1, -1], [5, -1], [5, 1], [1, 1]
c1 = hm.create_contour([p1, p2, p3, p4, p5, p1], [0, 1, 2, 3, 4])
g1 = hm.add_triangle_grid(p1, p2, p5, 5)
g2 = hm.add_unf_rect_grid(p2, p4, 10, 5)
g = hm.unite_grids(g1, [(g2, 0)])
g = hm.exclude_contours(g, c1, 'outer')
hm.export_grid_gmsh(g, "g1.msh")
hm.add_boundary_type(1, "b1")
hm.add_boundary_type(3, "b3")
g2 = hm.import_grid_gmsh("g1.msh")
hm.export_grid_gmsh(g2, "g2.msh")
hmdbg.check_ascii_file(8108300542786763763, "g1.msh")
hmdbg.check_ascii_file(17108615107233824167, "g2.msh")
hm.remove_all()

print "import from fluent msh file"
g1 = hm.import_grid_msh("../external_files/input1.msh")
hmdbg.check(hm.info_grid(g1) ==
            {'cell_types': {3: 182, 4: 1004},
             'Nnodes': 1196, 'Nedges': 2382, 'Ncells': 1186})
hmdbg.check(len(hm.registered_btypes()) == 6)
hm.export_grid_msh(g1, "g1.msh")
hmdbg.check_ascii_file(612863133040571191, "g1.msh")
hm.remove_all()

print "export 3d to fluent"
tribnd = hm.add_boundary_type(1, "triangle")
circbnd = hm.add_boundary_type(2, "circle")
botbnd = hm.add_boundary_type(3, "bottomz")
topbnd = hm.add_boundary_type(4, "topz")
g1 = hm.add_triangle_grid([0, 0], [1, 0], [0.5, 0.5], 10)
g2 = hm.add_unf_circ_grid([0.5, 1], 1.0, 30, 10, 1, False)
hm.set_boundary_type(g1, tribnd)
hm.set_boundary_type(g2, circbnd)
g3 = hm.unite_grids(g2, [(g1, 0.05)])
g4 = hm.extrude_grid(g3, [0, 0.05, 0.2, 0.5, 0.65, 0.7], botbnd, topbnd)
hm.export3d_grid_vtk(g4, None, "c1.vtk")
hmdbg.check_ascii_file(7889578680359330313, "c1.vtk", 'dev')
hm.export3d_grid_msh(g4, "g1.msh")
hmdbg.check_ascii_file(450400077272399620, "g1.msh", 'dev')

print "export 3d to fluent with periodic condition"
hm.export3d_grid_msh(g4, "g1.msh", [botbnd, topbnd, [0, 0, 0.0], [0, 0, 0.7]])
hmdbg.check_ascii_file(12326414976139824718, "g1.msh", 'dev')
hm.remove_all()

print "export contour to tecplot"
hm.add_boundary_type(1, "bbot")
hm.add_boundary_type(7, "bright")
c1 = hm.add_rect_contour([0, 0], [1, 1], [1, 0, 3, 7])
hm.export_contour_tecplot(c1, "c1.dat")
hmdbg.check_ascii_file(10478281133073733973, "c1.dat")

print "export 2d grid to tecplot"
g1 = hm.add_unf_circ_grid([1, 2], 10, 20, 8, 0.5, False)
g2 = hm.add_unf_ring_grid([-8, -7], 2, 5, 10, 3, 1.3)
hm.set_boundary_type(g1, 2)
hm.set_boundary_type(g2, 3)
g3 = hm.unite_grids(g1, [(g2, 0.5)], True, zero_angle_approx=180)
hm.export_contour_tecplot(g3, "c1.dat")
hmdbg.check_ascii_file(13195362596084466346, "c1.dat", 'dev')
hm.export_grid_tecplot(g3, "c1.dat")
hmdbg.check_ascii_file(7412218476507145895, "c1.dat", 'dev')

print "export 3d grid to tecplot"
g4 = hm.extrude_grid(g3, [0, 1, 2, 3, 5, 8, 9], 0, 10)
hm.export3d_grid_tecplot(g4, "c1.dat")
hmdbg.check_ascii_file(16956269303881327848, "c1.dat", 'dev')

g1 = hm.add_unf_rect_grid([0, 0], [10, 10], 10, 10)
c1 = hm.create_contour([[-2, 4], [8, -2], [4, 12], [-2, 4]])
g1 = hm.exclude_contours(g1, c1, "inner")
g2 = hm.extrude_grid(g1, [0, 0.2, 0.5, 0.9, 1.0], 1, 2, 3)
hm.export3d_grid_tecplot(g2, "c2.dat")
hmdbg.check_ascii_file(9823534259060690640, "c2.dat")

print "export with additional fields"
hm.remove_all()
g1 = hm.add_unf_rect_grid([0, 0], [10, 10], 10, 10)
c1 = hm.create_contour([[-1, -1], [12, 12], [12, -1], [-1, -1]])
g1 = hm.exclude_contours(g1, c1, "inner")
g2 = hm.extrude_grid(g1, [0, 0.2, 0.5, 0.9, 1.0], 1, 2, 3)
hm.export_grid_hmg(g1, "g1.hmg", afields=['cell-vertices',
                                          'cell-edges'])
hm.export3d_grid_hmg(g2, "g2.hmg", afields=['cell-vertices',
                                            'cell-faces',
                                            'face-vertices',
                                            'linfem'])
hmdbg.check_ascii_file(8700456673301720132, "g1.hmg")
hmdbg.check_ascii_file(5039114778377738832, "g2.hmg")
