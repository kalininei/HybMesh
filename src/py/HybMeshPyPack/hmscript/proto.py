from HybMeshPyPack import basic, com
import HybMeshPyPack.com.gridcom
from HybMeshPyPack.basic.geom import Point2
from HybMeshPyPack.hmscript import flow


# Prototype grids
def AddUnfRectGrid(p0, p1, nx, ny):
    """ add rectangular grid
        p0, p1 -- corner points as [x0, y0], [x1, y1]
        nx, ny -- partition in x, y direction
        returns grid identifier
    """
    c = com.gridcom.AddUnfRectGrid({"p0": Point2(*p0),
        "p1": Point2(*p1), "nx": nx, "ny": ny})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def AddUnfCircGrid(p0, rad, na, nr, coef=1, is_trian=True):
    """ add circular grid
        p0 - center coordinate as [x, y]
        rad - radius
        na, nr - partitions in arc and radius directions
        coef - refinement coefficient
            coef = 1: equidistant radius division
            coef > 1: refinement towards center of circle
            0< coef <1:  refinement towards outer arc
        is_trian - True if center should be triangulated
    """
    c = com.gridcom.AddUnfCircGrid({
        "p0": Point2(*p0),
        "rad": rad,
        "na": na, "nr": nr,
        "coef": coef,
        "is_trian": is_trian})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def AddUnfRingGrid(p0,
                radinner, radouter,
                na, nr, coef=1.0):
    """ add ring grid
        p0 - center coordinates as [x, y]
        radinner, radouter - inner and outer radii
        na, nr - arc and radius partition
        coef - refinement coefficient
            coef = 1: equidistant radius division
            coef > 1: refinement towards inner arc
            0< coef <1:  refinement towards outer arc
    """
    c = com.gridcom.AddUnfRingGrid({
        "p0": Point2(*p0),
        "radinner": radinner, "radouter": radouter,
        "na": na, "nr": nr,
        "coef": coef})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


# Contours
