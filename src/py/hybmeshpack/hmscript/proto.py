from hybmeshpack import com
from hybmeshpack.basic.geom import Point2, angle_3pnt
from hybmeshpack.hmscript import flow
from hybmeshpack.hmscript import ExecError
import math


# Prototype grids
def add_unf_rect_grid(p0, p1, nx, ny):
    """Builds rectangular grid

    Args:
       p0, p1 (list-of-float): bottom left, top right points as [x, y] list

       nx, ny (int): partition in x, y direction

    Returns:
       created grid identifier

    """
    c = com.gridcom.AddUnfRectGrid({"p0": Point2(*p0),
                                    "p1": Point2(*p1),
                                    "nx": nx,
                                    "ny": ny})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def add_unf_circ_grid(p0, rad, na, nr, coef=1, is_trian=True):
    """Builds circular grid

    Args:
       p0 (list-of-float): center coordinate as [x, y]

       rad (float): radius

       na, nr (int): partitions of arc and radius respectively

    Kwargs:
       coef (float): refinement coefficient::

         coef = 1: equidistant radius division
         coef > 1: refinement towards center of circle
         0 < coef < 1: refinement towards outer arc

       is_trian (bool): True if center cell should be triangulated

    Returns:
       created grid identifier

    """
    c = com.gridcom.AddUnfCircGrid({
        "p0": Point2(*p0),
        "rad": rad,
        "na": na, "nr": nr,
        "coef": coef,
        "is_trian": is_trian})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def add_unf_ring_grid(p0, radinner, radouter,
                      na, nr, coef=1.0):
    """Builds ring grid

    Args:
       p0 (list-of-float): center coordinates as [x, y]

       radinner, radouter (float): inner and outer radii

       na, nr (int): arc and radius partition respectively

    Kwargs:
       coef (float): refinement coefficient::

         coef = 1: equidistant radius division
         coef > 1: refinement towards center of circle
         0 < coef < 1: refinement towards outer arc

    Returns:
       created grid identifier

    """
    c = com.gridcom.AddUnfRingGrid({
        "p0": Point2(*p0),
        "radinner": radinner, "radouter": radouter,
        "na": na, "nr": nr,
        "coef": coef})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def add_triangle_grid(p0, p1, p2, nedge):
    """Creates structured grid in triangle area

    Args:
       ``p0, p1, p2``: triangle vertices in [x, y] format

       ``nedge`` (int): partition of triangle edges

    Returns:
       identifier of newly created grid

    Resulting grid will contain quadrangle cells everywhere except
    area near ``p0``-``p2`` edge where triangle cells will be built.
    """
    v = [Point2(*p0), Point2(*p1), Point2(*p2)]
    c = com.gridcom.AddTriGrid({
        "vertices": v, "nedge": nedge})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def add_custom_rect_grid(algo, left, bottom, right=None, top=None):
    """ Creates rectangular grid on the basis of four lines using
    contour vertices for partitioning.

    :param str algo: Algorithms of building:

       * ``'linear'`` - connects respective points of opposite
         contours by straight lines. If right/top contours are defined
         they should have same number of vertices as left/bottom.

    :param left:
    :param bottom:
    :param right:
    :param top: identifiers of base line segments.
       **right** and **top** could be ``None``. If so right and top
       boundaries will be created by translation of **left** and **bottom**.

    :return: new grid identifier
    :raise: hmscript.ExecError, ValueError

    If given contours are not properly connected then program will try
    to connect it with priority order:

    1) set left contour
    2) add bottom contour to left contour bottom point
    3) add top contour to left contour top point
    4) add right contour to bottom contour right point
    5) stretch right contour so its upper point fits top contour right point.

    """
    # call
    args = {'algo': algo,
            'left': left,
            'right': right,
            'bot': bottom,
            'top': top}
    c = com.gridcom.AddCustomRectGrid(args)
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('custom rectangular grid')


# Contour prototypes
def add_rect_contour(p0, p1, bnd=0):
    """Adds four point closed rectangular contour

    Args:
       p0, p1 (list-of-floats): bottom left and top right coordinates of
       the contour

    Kwargs:
       bnd (boundary identifier): single or list of 4 boundary
       identifiers (bottom, left, top, right) for contour segments.
       With the default value no boundary types will be set.

    Returns:
       Contour identifier

    """
    if isinstance(bnd, list):
        b = bnd[0:4]
    else:
        b = [bnd, bnd, bnd, bnd]

    c = com.contcom.AddRectCont({"p0": Point2(*p0),
                                 "p1": Point2(*p1),
                                 "bnds": b})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def add_circ_contour(p0, rad, n_arc, bnd=0):
    """Adds circle contour from given center and radius

    Args:
       p0 (list-of-floats): circle center

       rad (float): circle radius

       n_arc (int): partition of circle arc

    Kwargs:
       bnd (boundary identifier): boundary identifier for contour.
       With the default value no boundary types will be set.

    Returns:
       Contour identifier

    """
    c = com.contcom.AddCircCont({"p0": Point2(*p0), "rad": rad,
                                 "na": n_arc, "bnd": bnd})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def add_circ_contour2(p0, p1, p2, n_arc, bnd=0):
    """Adds circle contour from given arc points

    Args:
       p0, p1, p2: circle arc points as [x, y]

       n_arc (int): partition of circle arc

    Kwargs:
       bnd (boundary identifier): boundary identifier for contour.
       With the default value no boundary types will be set.

    Returns:
       Contour identifier

    """
    p0, p1, p2 = map(float, p0), map(float, p1), map(float, p2)
    try:
        xb, yb = p1[0] - p0[0], p1[1] - p0[1]
        xc, yc = p2[0] - p0[0], p2[1] - p0[1]
        A11, A12 = 2.0 * xb, 2.0 * yb
        A21, A22 = 2.0 * xc, 2.0 * yc
        B1, B2 = xb * xb + yb * yb, xc * xc + yc * yc
        d = A11 * A22 - A12 * A21
        I11, I12, I21, I22 = A22 / d, -A12 / d, -A21 / d, A11 / d
        cx = I11 * B1 + I12 * B2
        cy = I21 * B1 + I22 * B2
        rad = math.sqrt((cx - xb) * (cx - xb) + (cy - yb) * (cy - yb))
    except:
        raise ValueError("Failed to build a circle with given parameters")
    return add_circ_contour([cx + p0[0], cy + p0[1]], rad, n_arc, bnd)


def add_circ_contour3(p0, p1, curv, n_arc, bnd=0):
    """Adds circle contour from given arc points and curvature

    Args:
       p0, p1: circle arc points

       curv (float): circle curvature. Equals ``1.0/radius``.

       n_arc (int): partition of circle arc

    Kwargs:
       bnd (boundary identifier): boundary identifier for contour.
       With the default value no boundary types will be set.

    Returns:
       Contour identifier

    In the resulting circle ``p0``-``p1`` arc
    with counterclockwise direction will be shorter then
    ``p1``-``p0`` arc.
    """
    p0, p1, curv = map(float, p0), map(float, p1), float(curv)
    try:
        xa, ya = p1[0] - p0[0], p1[1] - p0[1]
        r = abs(1.0 / curv)
        a, b, c = -2.0 * xa, -2.0 * ya, xa * xa + ya * ya
        s = a * a + b * b
        x0, y0 = -a * c / s, -b * c / s
        d = r * r - c * c / s
        mult = math.sqrt(d / s)
        cx1 = x0 + b * mult
        cy1 = y0 - a * mult
        cx2 = x0 - b * mult
        cy2 = y0 + a * mult
        a1 = angle_3pnt(Point2(0.0, 0.0), Point2(cx1, cx2), Point2(xa, ya))
        if (a1 < math.pi):
            cx, cy = cx1, cy1
        else:
            cx, cy = cx2, cy2
    except:
        raise ValueError("Failed to build a circle with given parameters")
    return add_circ_contour([cx + p0[0], cy + p0[1]], r, n_arc, bnd)
