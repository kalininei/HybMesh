from hybmeshpack import hmscript as hm

hm.check_compatibility("0.3.0", 2)
ncirc = 256

# body
ntri, nspan = 12, 60
g3l = hm.add_triangle_grid([0, 0], [-0.5, 0.5], [0, 1], ntri)
g3r = hm.add_triangle_grid([5, 0], [5.5, 0.5], [5, 1], ntri)
g4 = hm.add_unf_rect_grid([0, 0], [5, 1], nspan, ntri)
gb = hm.unite_grids(g4, [(g3l, 0), (g3r, 0)])
circ1 = hm.add_circ_contour3([6, 0], [-6, 0], 0.08, ncirc)
circ2 = hm.add_circ_contour3([-6, 0], [6, 0], 0.08, ncirc)
bodyc = hm.clip_domain(circ1, circ2, "intersection")
body = hm.map_grid(gb, bodyc, [[-0.5, 0.5], [5.5, 0.5]],
                   [[-6, 0], [6, 0]], snap="no")

# tail
circ1 = hm.add_circ_contour3([0, 0], [-2, 2], 0.2, ncirc)
circ2 = hm.add_circ_contour3([-2.5, -2.3], [0, 0], 0.1, ncirc)
circ3 = hm.add_circ_contour3([-1.5, -1.8], [-2, 2], 0.2, ncirc)
d1 = hm.clip_domain(circ1, circ2, "intersection")
tailc = hm.clip_domain(d1, circ3, "difference")
g3 = hm.add_triangle_grid([-1, 0], [0, 2], [1, 0], 15)
tail = hm.map_grid(
    g3, tailc,
    [[-1, 0], [0, 2], [1, 0]],
    [[-2, 2], [0, 0], [-2.5, -2.3]])
hm.move_geom([tail], -5.7, 0)

#upper fin
circ1 = hm.add_circ_contour3([0.5, 1], [-0.7, 3.0], 0.4, ncirc)
circ2 = hm.add_circ_contour3([-1.5, 1], [-0.7, 3.0], 0.3, ncirc)
d1 = hm.clip_domain(circ1, circ2, "difference")
fin1c = hm.clip_domain(d1, body, "difference")
fin1 = hm.map_grid(
    g3, fin1c,
    [[1, 0], [0, 2], [-1, 0]],
    [[-0.7, 3], [-1.5, 1.3], [0.5, 1.3]], snap="add_vertices")

#side fin
fin2 = hm.add_triangle_grid([2.7, -0.6], [1.2, -1.0], [1.5, -2.8], 15)

#back fin
circ1 = hm.add_circ_contour2([-5.1, -1], [-5.0, 0], [-5.1, 1], ncirc)
circ2 = hm.add_circ_contour2([-5.1, -1], [-4, 0], [-5.1, 1], ncirc)
d1 = hm.clip_domain(circ2, circ1, "difference")
fin34c = hm.clip_domain(d1, body, "difference")
g3 = hm.add_triangle_grid([1, 0], [0, 0], [0, 1], 5)
fin3 = hm.map_grid(
    g3, fin34c,
    [[1, 0], [0, 0], [0, 1]],
    [[-4, 0.5], [-5.1, 0.3], [-5.1, 1]])
fin4 = hm.map_grid(
    g3, fin34c,
    [[0, 0], [1, 0], [0, 1]],
    [[-5.2, -0.3], [-5.1, -1], [-4.1, -0.7]])

#jaws
jcont = hm.create_contour([[4, -0.3], [8, -1.3], [8, -0.3], [4, -0.3]])
body = hm.exclude_contours(body, jcont, "inner")
btopt = hm.BoundaryGridOptions(
    body, [0, 0.07], "left", 0.1,
    start_point=[5, -0.5], end_point=[5.45, -0.36], project_to="corner")
btopt.range_angles[0] = 0
jgrid = hm.build_boundary_grid(btopt)

#gills
circ1 = hm.add_circ_contour3([2.7, -0.4], [2.7, 0.4], 0.4, ncirc)
circ2 = hm.add_circ_contour3([2.7, -0.4], [2.7, 0.4], 1.2, ncirc)
gillc = hm.clip_domain(circ2, circ1, "difference")
btopt = hm.BoundaryGridOptions(gillc, [0, 0.1], "right",
                               bnd_step=0.1, bnd_stepping="const")
gills = hm.build_boundary_grid(btopt)
[gills2] = hm.copy_geom(gills)
[gills3] = hm.copy_geom(gills)
hm.move_geom(gills2, -0.5, 0)
hm.move_geom(gills3, 0.5, 0)
gills = hm.unite_grids(gills, [(gills2, 0), (gills3, 0)])

#eye
eye = hm.add_unf_ring_grid([4.7, 0.3], 0.05, 0.1, 5, 1)

#coupling upper fin
hm.heal_grid(fin1, 10)
btopt = hm.BoundaryGridOptions(
    fin1, [0, 0.1], "right", bnd_stepping="no",
    start_point=[-1.2, -1.4], end_point=[0.5, 1.5])
tmp = hm.build_boundary_grid([btopt])
fin1 = hm.unite_grids(fin1, [(tmp, 0.0)])
shark = hm.unite_grids(body, [(fin1, 0.1)], zero_angle_approx=10)

#coupling backward fins
hm.heal_grid([fin3, fin4], 10)
btopt = hm.BoundaryGridOptions(
    fin3, [0, 0.1], "right", bnd_stepping="no",
    start_point=[-5.1, 0.4], end_point=[-4.3, 0.7])
tmp1 = hm.build_boundary_grid([btopt])
btopt = hm.BoundaryGridOptions(
    fin4, [0, 0.1], "right", bnd_stepping="no",
    start_point=[-4.3, -0.7], end_point=[-5.1, -0.4])
tmp2 = hm.build_boundary_grid([btopt])
fin3 = hm.unite_grids(fin3, [(tmp1, 0.0)])
fin4 = hm.unite_grids(fin4, [(tmp2, 0.0)])
shark = hm.unite_grids(shark, [(fin3, 0.1), (fin4, 0.1)], zero_angle_approx=10)

#coupling all other
shark = hm.unite_grids(shark, [(gills, 0.1), (tail, 0.1),
                               (eye, 0.1), (jgrid, 0.1), (fin2, 0.15)],
                       empty_holes=True, zero_angle_approx=10)
hm.heal_grid(shark)
print hm.skewness(shark)

hm.export_grid_vtk(shark, "shark.vtk")
