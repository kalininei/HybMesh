from hybmeshpack import hmscript as hm
from hybmeshpack.hmscript._dbg import check_ascii_file, checkdict, check
hm.check_compatibility("0.4.5")

print "add_unf_rect_grid"
g1 = hm.add_unf_rect_grid([0.01, -0.2], [1, 1], 15, 10)
checkdict(hm.info_grid(g1), {'Ncells': 150})
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(6298399812325595490, "g1.vtk")

g1 = hm.add_unf_rect_grid([0, 0], [1, 1], 15, 10, custom_x=0.5)
checkdict(hm.info_grid(g1), {'Ncells': 20})
g1 = hm.add_unf_rect_grid([1, 1], [2, 2], 15, 10, custom_y=0.2)
checkdict(hm.info_grid(g1), {'Ncells': 75})
g1 = hm.add_unf_rect_grid([1, 1], [2, 2], 3, 3, custom_x=0.1, custom_y=0.2)
checkdict(hm.info_grid(g1), {'Ncells': 50})
g1 = hm.add_unf_rect_grid([1, 1], [2, 2], custom_x=[0, 1], custom_y=[0, 1])
checkdict(hm.info_grid(g1), {'Ncells': 1})
g1 = hm.add_unf_rect_grid([1, 1], [2, 2], custom_x=[1, 1.9, 2],
                          custom_y=[-1, -0.2, 2.2])
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(11564798493749504112, "g1.vtk")

print "add_unf_circ_grid"
g1 = hm.add_unf_circ_grid([1, 1], 1, 13, 3, is_trian=True)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(8234138603355039682, "g1.vtk")
g1 = hm.add_unf_circ_grid([0, 0], 1, 5, 5, 0.3, is_trian=False)
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(7770728784250610107, "g1.vtk")
g1 = hm.add_unf_circ_grid([0, 0], 1, 6, 6, custom_rads=0.5)
checkdict(hm.info_grid(g1), {'Ncells': 12})
g1 = hm.add_unf_circ_grid([0, 0], 1, 6, 6, custom_rads=[0.5, 0.9, 1.5])
checkdict(hm.info_grid(g1), {'Ncells': 18})

g1 = hm.add_unf_circ_grid([0, 0], 2, 6, 6, custom_archs=[0.1])
checkdict(hm.info_grid(g1), {'Ncells': 756})

g1 = hm.add_unf_circ_grid([0, 0], 2, 6, 6,
                          custom_archs=[0.1, 0.2, 0.4, 0.7, 1])
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(2989160927046432181, "g1.vtk")

g1 = hm.add_unf_circ_grid([0, 0], 2, 6, 6,
                          custom_rads=[0.1, 0.2, 0.3, 0.4],
                          custom_archs=[45, 90, 180, 270, 405])
hm.export_grid_vtk(g1, "g1.vtk")
check_ascii_file(1918507880403773791, "g1.vtk")

print "partition_segment"
p = hm.partition_segment(1, 2, 0.1, 0.5)
check(round(sum(p), 2) == 7.02)
p = hm.partition_segment(1, 2, 0.1, 0.1, [1.5, 0.4])
check(round(sum(p), 2) == 9.)

print "hexagonal meshes"
g1 = hm.add_unf_hex_grid([[0, 0], 1], 0.1)
checkdict(hm.info_grid(g1), {'Ncells': 127})
g1 = hm.add_unf_hex_grid([[0, 0], 1], 0.088)
checkdict(hm.info_grid(g1), {'Ncells': 169})
g1 = hm.add_unf_hex_grid([[0, 0], 1], 0.11)
checkdict(hm.info_grid(g1), {'Ncells': 127})
g1 = hm.add_unf_hex_grid([[0, 0], 1], 1)
checkdict(hm.info_grid(g1), {'Ncells': 7})
g1 = hm.add_unf_hex_grid([[0, 0], 1], 2)
checkdict(hm.info_grid(g1), {'Ncells': 7})
g1 = hm.add_unf_hex_grid([[0, 0], 1], 0.1, True)
checkdict(hm.info_grid(g1), {'Ncells': 127})
g1 = hm.add_unf_hex_grid([[0, 0], 1], 5, True)
checkdict(hm.info_grid(g1), {'Ncells': 7})

g1 = hm.add_unf_hex_grid([[0, 0], [2, 1]], 0.1)
checkdict(hm.info_grid(g1), {'Ncells': 98})
g1 = hm.add_unf_hex_grid([[0, 0], [2, 1]], 0.09)
checkdict(hm.info_grid(g1), {'Ncells': 112})
g1 = hm.add_unf_hex_grid([[0, 0], [2, 1]], 1)
checkdict(hm.info_grid(g1), {'Ncells': 5})
g1 = hm.add_unf_hex_grid([[0, 0], [2, 1]], 3)
checkdict(hm.info_grid(g1), {'Ncells': 2})
g1 = hm.add_unf_hex_grid([[0, 0], [2, 1]], 0.1, True)
checkdict(hm.info_grid(g1), {'Ncells': 91})
g1 = hm.add_unf_hex_grid([[0, 0], [2, 1]], 1, True)
checkdict(hm.info_grid(g1), {'Ncells': 2})
g1 = hm.add_unf_hex_grid([[0, 0], [2, 1]], 5, True)
checkdict(hm.info_grid(g1), {'Ncells': 2})

