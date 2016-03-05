from hybmeshpack import hmscript
import copy

# START OF EXAMPLE
# create source contour for a boundary grid
cont = hmscript.add_rect_contour([0, 0], [1, 1])
# basic options for horizontal segments
# using default direction="left" to build grid inside the source contour
opvert = hmscript.BoundaryGridOptions(
    cont,
    partition=[0, 0.01, 0.02, 0.03],
    bnd_step=0.06)
# basic options for vertical segments
# It differs from horizontal segments by cross and lengthwise
# boundary grid partition
ophoriz = hmscript.BoundaryGridOptions(
    cont,
    partition=[0, 0.01, 0.02, 0.03, 0.05, 0.1],
    bnd_step=0.02)
# option for bottom segment. Using deepcopy so that changes in one
# inctance do not affect the others.
op1 = copy.deepcopy(ophoriz)
op1.start_point = [0, 0]
op1.end_point = [1, 0]
# option for left segment
op2 = copy.deepcopy(opvert)
op2.start_point = [1, 0]
op2.end_point = [1, 1]
# option for top segment
op3 = copy.deepcopy(ophoriz)
op3.start_point = [1, 1]
op3.end_point = [0, 1]
# option for right segment
op4 = copy.deepcopy(opvert)
op4.start_point = [0, 1]
op4.end_point = [0, 0]
# building boundary grid
bgrid = hmscript.build_boundary_grid([op1, op2, op3, op4])

# END OF EXAMPLE


print "build_boundary_grid example"
if abs(hmscript.domain_area(bgrid) - 0.248)>1e-7:
    raise Exception
