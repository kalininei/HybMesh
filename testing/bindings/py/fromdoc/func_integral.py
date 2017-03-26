# ===== <code_preproc remove vvv>
import sys
# This line was added to find Hybmesh.py that is not located
# in the current directory or package subdirectory.
# Normally this should be omitted.
sys.path.append("../../../../build/bindings/py")
# ===== <code_preproc remove ^^^>
import math
from Hybmesh import Hybmesh

# target function: Gaussian hill with the center at [0, 0]
def expfun(x, y):
    return math.exp(-(x*x + y*y)/(0.25))


# ===== <code_preproc remove vvv>
# set paths
Hybmesh.hybmesh_exec_path = "../../../../src/py/"
Hybmesh.hybmesh_lib_path = "../../../../build/bin/"
# ===== <code_preproc remove ^^^>

# create Hybmesh connection under with block to guarantee
# its destruction after all
with Hybmesh() as hm:

    # loop over different outer step sizes
    for h in [0.2, 0.1, 0.05]:

        # create rectangle-in-cirlce grid prototype
        g = hm.add_circ_rect_grid([0, 0], 1.0, h)

        # get vertices as plain 1D array
        vertices = g.raw_vertices()

        # get cell-vertices table as plain array.
        # Here we already know that grid contains only quadrangular cells,
        # hence there is no need to call g.raw_tab("cell_dim") to know
        # how the plain cell_vert array should be subdivided.
        cell_vert = g.raw_tab("cell_vert")

        # calculating integral as the sum of cell areas multiplied by
        # cell center function values.
        result = 0
        it = iter(cell_vert)
        for ind in zip(it, it, it, it):
            # ind contains four vertex indices for current cell
            # x, y - are respective x and y coordinates
            x = list(map(lambda i: vertices[2*i], ind))
            y = list(map(lambda i: vertices[2*i+1], ind))
            # calculate function value at approximate cell center
            f = expfun(sum(x)/4.0, sum(y)/4.0)
            # calculate area
            x1, x2, x3 = x[1] - x[0], x[2] - x[0], x[3] - x[0]
            y1, y2, y3 = y[1] - y[0], y[2] - y[0], y[3] - y[0]
            area = 0.5*(x1*y2-y1*x2 + x2*y3-y2*x3)
            result += f * area

        print("h = {}: integral = {}".format(h, result))
