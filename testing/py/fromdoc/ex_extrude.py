from hybmeshpack import hmscript
from hybmeshpack.hmscript._dbg import check_ascii_file
import math
global math

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# register boundary types
global bx0, bx1, by0, by1, bz0, bz1, bcustom
bx0 = hmscript.add_boundary_type(1, "x0-boundary")
bx1 = hmscript.add_boundary_type(2, "x1-boundary")
by0 = hmscript.add_boundary_type(3, "y0-boundary")
by1 = hmscript.add_boundary_type(4, "y1-boundary")
bz0 = hmscript.add_boundary_type(5, "z0-boundary")
bz1 = hmscript.add_boundary_type(6, "z1-boundary")
bcustom = hmscript.add_boundary_type(7, "bcustom")

# unit square 2D grid
square = hmscript.add_unf_rect_grid([0, 0], [1, 1], 5, 5)


# assign boundary types to 2D geometry so it could be inherited by 3D object
def assign_boundary2d(x0, y0, x1, y1, b):
    if x0 == 0 and x1 == 0:
        return bx0
    if x0 == 1 and x1 == 1:
        return bx1
    if y0 == 0 and y1 == 0:
        return by0
    if y0 == 1 and y1 == 1:
        return by1

hmscript.set_boundary_type(square, bfun=assign_boundary2d)

# calculate z coordinates with sine refinement towards z=0
minz, maxz, nz = 0.0, 1.0, 10


# refinement function [0, 1] -> [0, 1]
def reffun(t):
    return 1.0 + math.sin(0.5 * math.pi * (t - 1))

zcoords = []
for i in range(nz + 1):
    t = reffun(float(i) / nz)
    zcoords.append(minz + (maxz - minz) * t)

# case 1: extrude without any boundary assignment
res1 = hmscript.extrude_grid(square, zcoords, 0, 0, 0)


# case 2: assign all boundary types:
#         x0, x1, y0, y1 are inherited from 2d geometry,
#         z0, z1 are defined explicitly
res2 = hmscript.extrude_grid(square, zcoords, bz0, bz1)


# case 3: same as case 2 but select central face within z=0
#         surface and assign it with bcustom boundary type
def assign_boundary3d_z0(x, y, b):
    if (x - 0.5)**2 + (y - 0.5)**2 < 1e-6:
        return bcustom
    else:
        return bz0

res3 = hmscript.extrude_grid(square, zcoords, assign_boundary3d_z0, bz1)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
print "extrude example"
hmscript.export3d_grid_vtk(res1, None, "res1.vtk")
hmscript.export3d_grid_vtk(res2, None, "res2.vtk")
hmscript.export3d_grid_vtk(res3, None, "res3.vtk")
check_ascii_file(15866305684735266182, "res1.vtk")
check_ascii_file(501441713051147881, "res2.vtk")
check_ascii_file(11110412153695859606, "res3.vtk")
