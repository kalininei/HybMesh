import numbers
import math
import copy
from hybmeshpack import com
from hybmeshpack.basic.geom import Point2, angle_3pnt
from hybmeshpack.hmscript import flow
from hybmeshpack.hmscript import ExecError


# Prototype grids
def add_unf_rect_grid(p0=[0, 0], p1=[1, 1], nx=3, ny=3,
                      custom_x=None, custom_y=None):
    """Builds rectangular grid.

    :param list-of-floats p0:

    :param list-of-floats p1: bottom left, top right points in [x, y] format.

    :param int nx:

    :param int ny: partition in x and y directions.

    :returns: created grid identifier

    """
    if custom_x is None:
        custom_x = []
    elif not isinstance(custom_x, list):
        custom_x = [custom_x]
    if custom_y is None:
        custom_y = []
    elif not isinstance(custom_y, list):
        custom_y = [custom_y]
    c = com.gridcom.AddUnfRectGrid({"p0": Point2(*p0),
                                    "p1": Point2(*p1),
                                    "nx": nx,
                                    "ny": ny,
                                    "custom_x": custom_x,
                                    "custom_y": custom_y})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


def add_unf_circ_grid(p0, rad=1.0, na=8, nr=4, coef=1.0, is_trian=True,
                      custom_rads=None, custom_archs=None):
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
    if custom_rads is None:
        custom_rads = []
    if custom_archs is None:
        custom_archs = []
    if not isinstance(custom_rads, list):
        custom_rads = [custom_rads]
    if not isinstance(custom_archs, list):
        custom_archs = [custom_archs]
    c = com.gridcom.AddUnfCircGrid({
        "p0": Point2(*p0),
        "rad": rad,
        "na": na, "nr": nr,
        "coef": coef,
        "is_trian": is_trian,
        "custom_r": custom_rads,
        "custom_a": custom_archs})
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


def add_unf_hex_grid(area, cell_radius, strict=False):
    """ Builds grid with regular hexagonal cells
    """
    simpar = [area[0][0], area[0][1]]
    if isinstance(area[1], list):
        simpar.append(area[1][0])
        simpar.append(area[1][1])
    else:
        simpar.append(area[1])
    args = {"area": simpar, "crad": cell_radius, "strict": strict}
    c = com.gridcom.AddUnfHexGrid(args)
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('build uniform hexagonal grid')


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


def stripe(cont, partition, tip='no', bnd=None):
    """ Build a structured grid to the both sides of contour line

    :param cont: closed or open contour identifier

    :param ascending-list-of-double partition: partition perpendicular
       to source contour

    :param str tip: stripe endings meshing algorithm

       * ``"no"`` - no grid at endings
       * ``"radial"`` - radial grid at endings

    :param float-or-list-of-floats bnd: boundary types for input grid.
       List of four values provides respective values for bottom, left,
       right, top sides of resulting grid with respect to contour direction.

    :return: grid identifier

    :raise: hmscript.ExecError, ValueError

    Horizontal partition is taken from contour partition.
    Vertical partition is given by user with ``partition`` list parameter.
    If it starts with non zero value then grid will not contain
    contour nodes as its vertices.

    """
    #checks
    if tip not in ['no', 'radial']:
        raise ValueError("Invalid tip option")
    if not isinstance(partition, list):
        raise ValueError("Bad vertical partition")
    for p in partition:
        if not isinstance(p, (int, float, long)) or p < 0:
            raise ValueError("Bad vertical partition")
    for i in range(1, len(partition)):
        if partition[i] <= partition[i - 1]:
            raise ValueError("Bad vertical partition")
    if len(partition) < 2:
        if len(partition) == 0 or partition[0] == 0:
            raise ValueError("Bad vertical partition")

    if bnd is None:
        bnd = [0]
    elif not isinstance(bnd, list):
        bnd = [bnd]
    if len(bnd) < 4:
        for i in range(len(bnd), 4):
            bnd.append(bnd[-1])
    bnd = [bnd[0], bnd[1], bnd[2], bnd[3]]
    for b in bnd:
        if not isinstance(b, int):
            raise ValueError("Bad bnd parameter")
    arg = {"source": cont, "partition": partition, "tip": tip, "bnd": bnd}
    c = com.gridcom.StripeGrid(arg)
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('stripe grid')


def _triquad(domain, constr, pts, tp):
    if not isinstance(domain, list):
        domain = [domain]
    if constr is not None:
        if not isinstance(constr, list):
            constr = [constr]
    else:
        constr = []
    if pts is None:
        p2 = []
    else:
        p2 = []
        for i in range(len(pts) / 2):
            p2.extend(pts[2 * i + 1])
            p2.append(pts[2 * i])
    arg = {"domain": domain, "constr": constr, "pts": p2}
    if (tp == '3'):
        c = com.gridcom.TriangulateArea(arg)
    elif (tp == '4'):
        c = com.gridcom.QuadrangulateArea(arg)
    elif (tp == 'pebi'):
        c = com.gridcom.PebiFill(arg)
    else:
        raise ValueError
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('unstructured fill')


def triangulate_area(domain, constr=None, pts=None):
    """Builds constrained triangulation within given domain

    :param domain: single or list of closed contours
        representing bounding domain

    :param constr: single or list of contours representing
        triangulation constraints

    :param pts: set of points in ``[len0, [x0, y0], ...]``
        format where ``x, y`` are coordinates of internal verticies
        which should be embedded into the resulting grid,
        ``len`` - size of adjacent cells

    :return: grid identifier

    Triangulation procedure takes ``domain`` and ``constr`` verticies as
    basic points to define grid cell sizes.
    Use :func:`partition_contour` procedure to adopt partition of these
    contours and control resulting cell sizes.

    It is guaranteed that no grid edge will cross ``constr`` contours.

    Crosses between constrain contours and domain are supported.
    """
    return _triquad(domain, constr, pts, '3')


def quadrangulate_area(domain, constr=None, pts=None):
    """Builds constrained quadrangulation within given domain

    :param domain: single or list of closed contours
        representing bounding domain

    :param constr: single or list of contours representing
        meshing constraints

    :param pts: set of points in ``[len0, [x0, y0], ...]``
        format where ``x, y`` are coordinates of internal verticies
        which should be embedded into the resulting grid,
        ``len`` - size of adjacent cells

    :return: grid identifier

    See :func:`triangulate_area` for parameters details.

    .. note::

       This procedure doesn't guarantee that all cells in the resulting
       grid will be quadrangular. If input geometry is not good enough
       some triangle cells may occur.

    """
    return _triquad(domain, constr, pts, '4')


def pebi_fill(domain, constr=None, pts=None):
    """Builds perpendicular bisector cells in given domain

    :param domain: single or list of closed contours
        representing bounding domain

    :param constr: single or list of contours representing
        meshing constraints

    :param pts: set of points in ``[len0, [x0, y0], ...]``
        format where ``x, y`` are coordinates of internal verticies
        which should be embedded into the resulting grid,
        ``len`` - size of adjacent cells

    :return: grid identifier

    Domain points as well as points of ``pts`` list and constrain contour
    nodes will be treated as pebi cell centers.
    Therefore boundary representation (but not domain area) of resulting grid
    will be different from passed  as ``domain`` parameter.

    Routine can produce concave cells (f.e. as a result of bad size
    control or near the concave domain boundary vertices).
    Use :func:`heal_grid` routine with ``convex_cells`` option to fix this.

    .. note::

       After the grid is build some optimisation procedures will be executed
       in order to get rid of short edges and possible self intersections.
       So the resulting grid will not be strictly of pebi type.
    """
    return _triquad(domain, constr, pts, 'pebi')


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
