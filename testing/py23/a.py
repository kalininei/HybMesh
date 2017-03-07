import sys
sys.path.append("../../build/bindings/py")
from Hybmesh import Hybmesh  # NOQA


def check_cond(cond):
    if (not cond):
        raise Exception("check failed")


def check_dims(v1, v2):
        if len(v1) != len(v2):
            raise Exception("check failed")
        for i in range(len(v2)):
            if v1[i] != v2[i]:
                raise Exception("check failed")


def good_callback(s1, s2, p1, p2):
    print("{} {} -- {} {}".format(s1, p1, s2, p2))
    return 0


def bad_callback(s1, s2, p1, p2):
    print("{} {} -- {} {}".format(s1, p1, s2, p2))
    return 0 if p1 < 0.7 else 1


Hybmesh.hybmesh_exec_path = "../../src/py/hybmesh.py"
Hybmesh.hybmesh_lib_path = "../../build/bin/"

with Hybmesh() as hm:
    g1 = hm.add_unf_rect_grid1([1, 2, 3], [2, 3, 4])
    dims = g1.dims()
    check_dims(dims, [9, 12, 4])

    c1 = hm.grid_bnd_to_contour(g1, True)
    c2 = hm.grid_bnd_to_contour(g1, False)
    check_dims(c1.dims(), [4, 4])
    check_dims(c2.dims(), [8, 8])
    pv = [None]*4
    pv[0] = [0, 0]
    pv[1] = [1, 2]
    pv[2] = [3, 1]
    pv[3] = [2, 0]
    c3 = hm.create_contour(pv, [0, 0, 1])
    c4 = hm.create_contour(pv, [3, 2, 1])
    check_dims(c3.dims(), [4, 3])
    check_dims(c4.dims(), [4, 3])
    c5 = hm.create_spline_contour(pv, [0, 0, 1])
    c6 = hm.create_spline_contour(pv, [3, 2, 1])
    check_dims(c5.dims(), [101, 100])
    check_dims(c6.dims(), [101, 100])
    c7 = hm.extract_subcontours(c6, pv)
    check_cond(len(c7) == 3)
    c8 = hm.add_rect_contour([1, 1], [10, 10], [0, 1, 2, 3])
    check_dims(c8.dims(), [4, 4])
    c9 = hm.add_circ_contour([5., 5.], 10, 32)
    check_dims(c9.dims(), [32, 32])
    ar9 = c9.domain_area()
    check_cond(abs(ar9 - 3.1415926*100) < 5)
    c10 = hm.simplify_contour(c2)
    check_dims(c10.dims(), [4, 4])
    c11 = hm.simplify_contour(c9, 91)
    check_dims(c11.dims(), [4, 4])
    c12 = hm.unite_contours([c8, c9])
    check_dims(c12.dims(), [36, 36])
    c13 = hm.decompose_contour(c12)
    check_cond(len(c13) == 2)
    c14 = hm.clip_domain(c13[0], c13[1], "difference")
    c15 = hm.add_rect_contour([0, 0], [1, 1])
    c16 = hm.add_rect_contour([10, 10], [11, 11])
    c17 = hm.clip_domain(c15, c16, "intersection")
    c18 = hm.clip_domain(c13[1], c13[1], "difference")
    check_dims(c17.dims(), [0, 0])
    check_dims(c18.dims(), [0, 0])
    c19 = hm.partition_contour_const(c15, 0.01)
    check_dims(c19.dims(), [400, 400])
    c20 = hm.partition_contour_const(c15, 0.01)
    check_dims(c19.dims(), [400, 400])
    c21 = hm.partition_contour_ref_points(
        c15, [0.01, 0.5], [[0, 0], [1, 1]],
        30, True, -1, [], [], [0, 0], [1, 1])
    check_dims(c21.dims(), [18, 18])
    c22 = hm.partition_contour_ref_lengths(
        c15, [0.01, 0.5], [0, 2],
        30, True, -1, [], [], [0, 0], [1, 1])
    check_dims(c21.dims(), c22.dims())
    c23 = hm.matched_partition(
        c15, 0.1, 0.5, [], [0.01], [[0.1, 0.1]])
    check_dims(c23.dims(), [74, 74])
    c24 = hm.partition_segment(0.5, 1.5, 0.1, 0.5, [0.01, 1.0])
    check_cond(len(c24) == 5 and c24[0] == 0.5)

    c25 = hm.create_contour([[0, 0], [1, 0], [2, 0]])
    c26 = hm.create_contour([[2, 0.1], [3, 0.1], [4, 1]])
    c27 = hm.connect_subcontours([c25, c26], [])
    c28 = hm.connect_subcontours([c25, c26], [], "yes")
    c29 = hm.connect_subcontours([c25, c26], [], "force")
    check_dims(c27.dims(), [5, 4])
    check_dims(c28.dims(), [4, 4])
    check_dims(c29.dims(), [5, 5])

    c30 = hm.add_unf_rect_grid1([1, 2, 3], [4, 5, 6, 7])
    check_dims(c30.dims(), [12, 17, 6])
    c31 = hm.add_unf_circ_grid([0, 0], 10, 32, 5, 1.2, False)
    check_dims(c31.dims(), [160, 288, 129])
    c32 = hm.add_unf_ring_grid([1, 1], 3, 7, 16, 2)
    check_dims(c32.dims(), [48, 80, 32])
    c33 = hm.add_unf_hex_grid_in_hex([-1, -1], 8, 1)
    c34 = hm.add_unf_hex_grid_in_rect([0, 0], [2, 3], 1)
    c35 = hm.add_triangle_grid([0, 0], [1, 0], [0, 1], 3, [1, 2, 3])
    check_dims(c33.dims(), [216, 306, 91])
    check_dims(c34.dims(), [28, 35, 8])
    check_dims(c35.dims(), [10, 15, 6])

    c36 = hm.create_contour([[0, 0], [0, 1]])
    c37 = hm.partition_contour_const(c36, 0.1)
    c38 = hm.create_contour([[0, 0], [2, -0.1]])
    c39 = hm.partition_contour_const(c38, 0.2)
    c40 = hm.add_custom_rect_grid("linear", c37, c39)
    c41 = hm.add_custom_rect_grid_htfi(
        c37, c39, None, None, [1, 1, 1, 0.8])
    hm.assign_callback(bad_callback)
    try:
        c42 = hm.add_custom_rect_grid("orthogonal", c37, c39)
    except Hybmesh.EUserInterrupt as e:
        print("User interrupt catched")
    hm.assign_callback(good_callback)
    c42 = hm.add_custom_rect_grid("orthogonal", c37, c39)
    hm.reset_callback()
    check_dims(c40.dims(), [121, 220, 100])
    check_dims(c40.dims(), c41.dims())
    check_dims(c40.dims(), c42.dims())
    try:
        c43 = hm.add_circ_rect_grid([0, 0], -1, 0.05)
    except Hybmesh.ERuntimeError as e:
        print("Runtime error catched")
        print(str(e))
    c43 = hm.add_circ_rect_grid([0, 0], 1, 0.05)
    c44 = hm.add_circ_rect_grid([0, 0], 1, 0.05, 1., 1., "orthogonal_rect")
    check_dims(c43.dims(), c44.dims())
    c45 = hm.stripe(c39, [0, 0.01, 0.02, 0.05])
    c46 = hm.triangulate_domain(c45, [], [], [], "3")
    c47 = hm.grid_bnd_to_contour(c45)
    c48 = hm.grid_bnd_to_contour(c46)
    check_cond(abs(c47.domain_area() - c48.domain_area()) < 1e-8)
    c49 = hm.add_rect_contour([0, 0], [1, 1])
    c50 = hm.partition_contour_const(c49, 0.1)
    c51 = hm.pebi_fill(c50, [], [0.5], [[0.5, 0.5]])
    c52 = hm.simple_boundary_grid(c51, [0, 0.01, 0.02])
    c53 = hm.simple_boundary_grid(c51, [0, 0.01, 0.02], "right")
    check_dims(c52.dims(), c53.dims())
    c54 = hm.exclude_contours(c53, [c51], "inner")
    c55 = hm.exclude_contours(c53, [c51], "outer")
    check_dims(c54.dims(), c53.dims())
    check_dims(c55.dims(), [0, 0, 0])
    c56x, c57x = [], []
    for i in range(11):
            c56x.append(float(i)/10.)
            c57x.append(0.3 + float(i)/30.)
    c56 = hm.add_unf_rect_grid1(c56x, c56x)
    c57 = hm.add_unf_rect_grid1(c57x, c57x)
    c58 = hm.unite_grids1(c56, c57, 0.1)
    c59 = hm.map_grid(c58, c57, [[0, 0]], [[0.3, 0.3]])
    hm.export_grid_vtk(c58, "c58.vtk");
    hm.heal_grid(c59, 30, 30)
    c60 = hm.exclude_contours(c58, [c57], "inner")
    c61 = hm.extrude_grid(c60, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
    c62 = hm.grid3_bnd_to_surface(c61)
    hm.assign_callback(good_callback)
    c63 = hm.tetrahedral_fill([c62])
    hm.reset_callback()
    c64 = hm.revolve_grid(c60, [0, 0], [0, 1],
                          [0, 45, 90, 180])
    c65 = hm.revolve_grid(c60, [0, 0], [0, 1],
                          [180, 270, 360])
    hm.stdout_verbosity(3)
    hm.export3d_grid_vtk(c64, "c64.vtk")
    hm.stdout_verbosity(0)
    hm.export3d_grid_vtk(c65, "c65.vtk")
    c66 = hm.merge_grids3(c64, c65)
    hm.export3d_grid_vtk(c66, "c66.vtk")
    c67 = c66.deepcopy()
    c67.scale(20, 20, 20, [0, 0, 0])
    c68 = hm.grid3_bnd_to_surface(c66)
    c69 = hm.grid3_bnd_to_surface(c67)
    check_cond(abs(c68.domain_volume()-125*c69.domain_volume()) < 1e-8)
    c60.rotate(20, [0, 0])
    c69.free()

    hm.add_boundary_type(1, "bnd1")
    hm.add_boundary_type(2, "bnd2")
    c70 = hm.add_unf_rect_grid([0, 0], [1, 1], 10, 10, [0, 1, 2, 3])
    hm.export_grid_msh(c70, "out.msh", [0], [2], [False])

    c71 = c70.skewness(-1)
    c72 = c70.skewness(0.1)
    check_cond(len(c71) == 100 and len(c72) == 0)

    c70.set_btypes_all(22)
    hm.export_contour_vtk(c70, "out.vtk")
    c73 = c70.raw_vertices()
    check_cond(len(c73) == 242)
    c74 = c70.raw_tab("bnd_bt")
    check_cond(len(c74) == 80)
    c75 = c70.raw_tab("bnd")
    check_cond(len(c75) == 40)
    c76 = c70.raw_tab("cell_edge")
    check_cond(len(c76) == 400)
    c77 = c70.raw_tab("edge_cell")
    check_cond(len(c77) == 440)
    c78 = c70.raw_tab("bt")
    check_cond(len(c78) == 220 and all(map(lambda x: x in [0, 22], c78)))
    c79 = c70.raw_tab("cell_dim")
    check_cond(len(c79) == 100 and all(map(lambda x: x == 4, c79)))
    c80 = c70.raw_tab("cell_vert")
    check_cond(len(c80) == 400)
    c70.set_btypes(1, c75[:3])
    hm.export_grid_tecplot(c70, "out.dat")

    c81 = hm.extrude_grid(c70, [0, 1, 2])
    c82 = hm.grid3_bnd_to_surface(c81)
    c83 = c81.raw_vertices()
    c84 = c82.raw_vertices()
    check_cond(len(c83) == 121*3*3)
    check_cond(len(c84) == (121*2 + 40)*3)
    c85 = c81.raw_tab("cell_vert")
    check_cond(len(c85) == 200*8)
    c86 = c82.raw_tab("face_vert")
    check_cond(len(c86) == (200+80)*4)

    c87 = hm.grid_bnd_to_contour(c70)
    c88 = c87.raw_vertices()
    check_cond(len(c88) == 12)
    c89 = c87.raw_tab("bt")
    check_cond(len(c89) == 6)

    c90 = hm.add_rect_contour([1, 1.2], [2, 2], [0, 1, 2, 3])
    c91 = hm.add_rect_contour([0, 0], [1.5, 1.5], [0, 1, 2, 3])
    c92 = hm.unite_contours([c90, c91])
    c93 = hm.decompose_contour(c92)
    check_cond(len(c93) == 3)
    c94 = hm.pick_contour([-1, -1], c93)
    c94.get_point([-1, -1])
    c94.get_point(None, [0.5, -1])
    c94.get_point(None, None, [1.0, -1])

    hm.stdout_verbosity(3)
    hm.export_grid_hmg([c55], "c55.hmg")
    hm.export_all_hmd("out.hmd", "ascii")
    hm.export_all_hmd("out1.hmd", "bin")
    hm.remove_all()
    hm.import_grid_hmg("out.hmd", "Grid2D_1")
    hm.import3d_grid_hmg("out.hmd", "Grid3D_1")
    hm.remove_all()
