from HybMeshPyPack import hmscript as hm
# import copy
import math
global hm, check


def check(cond):
    if not cond:
        raise Exception("Test Failed")


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
        check(info['cells_types'][k] == ct[k])


def check_zero(a):
    check(abs(a) < 1e-8)

# print "simple bgrid around a square"
# c1 = hm.add_rect_cont([0, 0], [1, 1], 5)
# op1 = hm.BoundaryGridOption(c1, [0, 0.01, 0.02], "left", 0.1)
# g1 = hm.build_boundary_grid(op1)
# op2 = hm.BoundaryGridOption(c1, [0, 0.01, 0.02], "right", 0.1)
# g2 = hm.build_boundary_grid(op2)

# check_grid(g1, 144, 240, 96, {4: 96})
# check_cont(g1, 96, 96, [56, 40], {0: 40, 5: 56})
# check_grid(g2, 144, 240, 96, {4: 96})
# check_cont(g2, 96, 96, [56, 40], {0: 56, 5: 40})

# print "set of options for a square source"
# ophoriz = hm.BoundaryGridOption(c1, [0, 0.01, 0.02, 0.03],
#                                 "left", 0.05)
# opvert = hm.BoundaryGridOption(c1, [0, 0.01, 0.02, 0.03, 0.05, 0.1],
#                                "left", 0.03)
# op1 = copy.deepcopy(ophoriz)
# op1.start_point, op1.end_point = [0, 0], [1, 0]
# op2 = copy.deepcopy(opvert)
# op2.start_point, op2.end_point = [1, 0], [1, 1]
# op3 = copy.deepcopy(ophoriz)
# op3.start_point, op3.end_point = [1, 1], [0, 1]
# op4 = copy.deepcopy(opvert)
# op4.start_point, op4.end_point = [0, 1], [0, 0]
# g3 = hm.build_boundary_grid([op1, op2, op3, op4])
# op1.direction = op2.direction = op3.direction = op4.direction = "right"
# g4 = hm.build_boundary_grid([op1, op2, op3, op4])

# check_cont(g3, 220, 220, [126, 94], {0: 94, 5: 126})
# check_cont(g4, 244, 244, [106, 138], {0: 138, 5: 106})
# check_grid(g3, 576, 1042, 466, {4: 466})
# check_grid(g4, 632, 1142, 510, {4: 510})

print "increasing angle test"
start, end, diff_ac, diff = 20, 350.0, 2, 50.0
angle = start
while angle < end:
    print "angle =", angle
    a = angle / 360.0 * 2 * math.pi
    cont = hm.create_contour([[-1, 0], [0, 0],
                              [-math.cos(a), math.sin(a)]])
    op = hm.BoundaryGridOption(cont, [0, 0.005, 0.01, 0.017, 0.027, 0.04],
                               'left', 0.01)
    grid = hm.build_boundary_grid(op)
    hm.export_grid_vtk(grid, ''.join(["angle_", str(angle), ".vtk"]))
    hm.remove_geom([cont, grid])
    angle += diff_ac if angle < op.range_angles[0] else diff
    break
