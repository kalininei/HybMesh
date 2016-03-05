from hybmeshpack import hmscript

# START OF EXAMPLE
# mark b1, b2 as globals so binfo function can find them
global b1, b2
# register boundary types
b1 = hmscript.add_boundary_type(1, "vertical")
b2 = hmscript.add_boundary_type(2, "horizontal")
# create rectangular contours
cont1 = hmscript.add_rect_contour([0, 0], [1, 1])
# mark contour using `btps`
hmscript.set_boundary_type(cont1, btps=[b2, b1, b2, b1])
# same effect using `bfun`
cont2 = hmscript.add_rect_contour([0, 0], [1, 1])
def binfo(x0, y0, x1, y1, bold):
    return b1 if x0 == x1 else b2
hmscript.set_boundary_type(cont2, bfun=binfo)
# END OF EXAMPLE


print "set_boundary_type example"
hmscript.export_contour_vtk(cont1, "_c1.vtk")
hmscript.export_contour_vtk(cont2, "_c2.vtk")
f1, f2 = open("_c1.vtk", "r"), open("_c2.vtk", "r")
if hash(f1.read()) != hash(f2.read()):
    raise Exception
if (hmscript.info_contour(cont1)['btypes'] != {1: 2, 2: 2}):
    raise Exception
