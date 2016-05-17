import numbers
import math
import copy
from hybmeshpack import com
from hybmeshpack.basic.geom import Point2, angle_3pnt
from hybmeshpack.hmscript import flow
from hybmeshpack.hmscript import ExecError


# Prototype grids
def add_unf_rect_grid(p0, p1, nx, ny):
    """Builds rectangular grid.

    :param list-of-floats p0:

    :param list-of-floats p1: bottom left, top right points in [x, y] format.

    :param int nx:

    :param int ny: partition in x and y directions.

    :returns: created grid identifier

    """
    c = com.gridcom.AddUnfRectGrid({"p0": Point2(*p0),
                                    "p1": Point2(*p1),
                                    "nx": nx,
                                    "ny": ny})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def add_unf_circ_grid(p0, rad, na, nr, coef=1.0, is_trian=True):
    """Builds circular grid.

    :param list-of-floats p0: center coordinate as [x, y]

    :param float rad: radius

    :param int na:

    :param int nr: partitions of arc and radius respectively

    :param float coef: refinement coefficient:
         * ``coef = 1``: equidistant radius division
         * ``coef > 1``: refinement towards center of circle
         * ``0 < coef < 1``: refinement towards outer arc

    :param bool is_trian: True if center cell should be triangulated

    :returns: created grid identifier

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

    :param list-of-floats p0: center coordinates as [x, y]

    :param float radinner:

    :param float radouter: inner and outer radii

    :param int na:

    :param int nr: arc and radius partition respectively

    :param float coef: refinement coefficient:
       * ``coef = 1``: equidistant radius division
       * ``coef > 1``: refinement towards center of circle
       * ``0 < coef < 1``: refinement towards outer arc

    :return: created grid identifier

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

    :param list-of-floats p0:

    :param list-of-floats p1:

    :param list-of-floats p2: triangle vertices in [x, y] format

    :param int nedge: partition of triangle edges

    :return: identifier of newly created grid

    Resulting grid will contain quadrangle cells everywhere except
    area near ``p0``-``p2`` edge where triangle cells will be built.
    """
    v = [Point2(*p0), Point2(*p1), Point2(*p2)]
    c = com.gridcom.AddTriGrid({
        "vertices": v, "nedge": nedge})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def add_custom_rect_grid(algo, left, bottom, right=None, top=None,
                         hermite_tfi_w=[1.0, 1.0, 1.0, 1.0]):
    """ Creates rectangular grid on the basis of four curvilinear contours
    using contour vertices for partition.
    See details in :ref:`custom_rect_grid`.

    :param str algo: Algorithms of building:

       * ``'linear'`` - connects respective points of opposite
         contours by straight lines.
       * ``'linear_tfi'`` - linear transfinite interpolation
       * ``'hermite_tfi'`` - hermite transfinite interpolation
       * ``'inverse_laplace'`` -
       * ``'direct_laplace'`` - connects points using solution of
         laplace equation with Dirichlet boundary conditions;
       * ``'orthogonal'`` - builds orthogonal grid based on **left** and
         **bottom** partition.
         Partitions of **right** and **top** are ignored.

    :param left:

    :param bottom:

    :param right:

    :param top: identifiers of curvilinear domain sides.
       **right** and **top** could be ``None``. If so right and top
       boundaries will be created by translation of **left** and **bottom**.

    :param list-of-floats hermite_tfi_w:
       perpendicularity weights
       for **left**, **bottom**, **right**, **top** contours respectively
       for **algo** = ``'hermite_tfi'``

    :return: new grid identifier

    :raise: hmscript.ExecError, ValueError

    """
    if algo not in ['linear', 'inverse_laplace', 'direct_laplace',
                    'orthogonal', 'linear_tfi', 'hermite_tfi']:
        raise ValueError("Unknown algorithm %s" % str(algo))
    if not isinstance(hermite_tfi_w, list) or len(hermite_tfi_w) != 4:
        raise ValueError("invalid hermite_tfi perpendicularity weights")
    # call
    args = {'algo': algo,
            'left': left,
            'right': right,
            'bot': bottom,
            'top': top,
            'her_w': copy.deepcopy(hermite_tfi_w)}
    c = com.gridcom.AddCustomRectGrid(args)
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('custom rectangular grid')


def add_circ_rect_grid(p0, rad, step, sqrside=1.0, rcoef=1.0, algo="linear"):
    """ Creates quadrangular cell grid in a circular area.
    See details in :ref:`circrect_grid`.

    :param list-of-floats p0: center point of circle area in [x, y] format.

    :param positive-float rad: radius of circle area.

    :param positive-float step: approximate partition step of the outer
       boundary.

    :param positive-float sqrside: side of the inner square normalized by
       the circle radius. Values greater than 1.4 are not allowed.

    :param positive-float rcoef: radius direction refinement of
       the ring part of the grid.
       Values less then unity lead to refinement towards outer boundary.

    :param str algo: Algorithms of assembling the ring part of the grid.

       * ``'linear'`` - use weighted approach for each ray partition
       * ``'laplace'`` - use algebraic mapping for
         building each 45 degree sector.
       * ``'orthogonal_circ'`` - build orthogonal grid keeping
         uniform grid at outer circle
       * ``'orthogonal_rect'`` - build orthogonal grid keeping
         uniform grid at inner rectangle

    :return: new grid identifier

    :raise: hmscript.ExecError, ValueError

    """
    # check input data
    if (not isinstance(p0, list) or len(p0) != 2 or
            not isinstance(p0[0], numbers.Real) or
            not isinstance(p0[1], numbers.Real)):
        raise ValueError("Invalid center point")
    if (not isinstance(rad, numbers.Real) or rad <= 0):
        raise ValueError("Invalid radius")
    if (not isinstance(step, numbers.Real) or step <= 0):
        raise ValueError("Invalid step")
    if (not isinstance(sqrside, numbers.Real) or
            sqrside <= 0 or sqrside > 1.4):
        raise ValueError("Invalid sqrside")
    if (not isinstance(rcoef, numbers.Real) or rcoef <= 0):
        raise ValueError("Invalid rcoef")
    if algo not in ["laplace", "linear", "orthogonal_rect", "orthogonal_circ"]:
        raise ValueError("Unknown algorithm")
    # call
    args = {'algo': algo,
            'p0': Point2(*p0),
            'rad': rad,
            'step': step,
            'sqrside': sqrside,
            'rcoef': rcoef}
    c = com.gridcom.AddCirc4Grid(args)
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('circ_rect_grid')


# Contour prototypes
def add_rect_contour(p0, p1, bnd=0):
    """Adds four point closed rectangular contour

    :param list-of-floats p0:

    :param list-of-floats p1: bottom left and top right coordinates of
       the contour

    :param bnd: single or list of 4 boundary
       identifiers (bottom, left, top, right) for contour segments.
       With the default value no boundary types will be set.

    :return: Contour identifier
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

    :param list-of-floats p0: circle center in [x, y] format

    :param float rad: circle radius

    :param int n_arc: partition of circle arc

    :param bnd: boundary identifier for contour.
       With the default value no boundary types will be set.

    :return: Contour identifier
    """
    c = com.contcom.AddCircCont({"p0": Point2(*p0), "rad": rad,
                                 "na": n_arc, "bnd": bnd})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def add_circ_contour2(p0, p1, p2, n_arc, bnd=0):
    """Adds circle contour from given arc points

    :param list-of-floats p0:

    :param list-of-floats p1:

    :param list-of-floats p2: circle arc points as [x, y] format

    :param int n_arc: partition of circle arc

    :param bnd: boundary identifier for contour.
       With the default value no boundary types will be set.

    :return:  Contour identifier

    """
    try:
        p0, p1, p2 = map(float, p0), map(float, p1), map(float, p2)
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

    :param list-of-floats p0:

    :param list-of-floats p1: circle arc points in [x, y] format

    :param float curv: circle curvature. Equals ``1.0/radius``.

    :param int n_arc: partition of circle arc

    :param bnd: boundary identifier for contour.
       With the default value no boundary types will be set.

    :return: Contour identifier

    In the resulting circle ``p0``-``p1`` arc
    with counterclockwise direction will be shorter then
    ``p1``-``p0`` arc.
    """
    try:
        p0, p1, curv = map(float, p0), map(float, p1), float(curv)
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
