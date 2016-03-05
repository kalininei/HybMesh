from hybmeshpack import hmscript

# START OF EXAMPLE

# source grid: circle with diameter = 5
g1 = hmscript.add_unf_circ_grid([0, 0], 5, 16, 5)
# outer domain border: square with edge length = 8
c1 = hmscript.add_rect_contour([-4, -4], [4, 4])
# inner domain border: triangle within c1
c2 = hmscript.create_contour([[-2, -2], [2, -2], [0, 2], [-2, -2]])
# multiply connected domain
c3 = hmscript.unite_contours([c1, c2])
# exclusion
g2 = hmscript.exclude_contours(g1, c3, "outer")

# END OF EXAMPLE


print "exclude_contours example"
if hmscript.info_grid(g2)['Ncells'] != 66:
    raise Exception
