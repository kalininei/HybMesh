" 2d contours construction procedures"
import math
from hybmeshpack import com
from hybmeshpack.hmscript import flow, hmscriptfun
from datachecks import (icheck, List, Bool, Point2D, ListOr1,
                        ZType, Grid2D, UInt, ACont2D, OneOf, Float)


@hmscriptfun
def grid_bnd_to_contour(gid, simplify=True):
    """ Extracts grid boundary to user contour.

    :param gid: grid identifier

    :param bool simplify: if true deletes all non-significant points. Otherwise
       resulting contour will contain all boundary grid edges.

    :returns: contour identifier
    """
    icheck(0, Grid2D())
    icheck(1, Bool())

    c = com.contcom.GridBndToContour({"grid_name": gid,
                                      "simplify": simplify})
    flow.exec_command(c)
    return c.added_contours2()[0]


@hmscriptfun
def create_contour(pnts, bnds=0):
    """ Create singly connected contour from sequence of points.

    :param list-of-list-of-floats pnts: sequence of points.
       If coordinates of first and last points are equal
       then contour is considered closed.

    :param single-or-list-of-boundary-identifiers bnds: boundary type for
       each contour segment or single identifier for the whole contour.

    :returns: contour identifier

    Example:
       >>> hmscript.create_contour([[0, 0], [1, 0], [1, 1], [0, 0]],
                                  [b1, b2, b3])
    """
    icheck(0, List(Point2D(), minlen=2))
    icheck(1, ListOr1(ZType()))

    b = bnds if isinstance(bnds, list) else [bnds]
    c = com.contcom.CreateContour({"points": pnts,
                                   "bnds": b})
    flow.exec_command(c)
    return c.added_contours2()[0]


@hmscriptfun
def create_spline_contour(pnts, bnds=0, nedges=100):
    """ Creates singly connected contour as a parametric cubic spline.

    :param list-of-list-of-floats pnts: sequence of points.
         If coordinates of first and last points are equal
         then resulting contour will be closed.

    :param single-or-list-of-boundary-identifiers bnds: boundary type for
         each contour segment bounded by **pnts**
         or single identifier for the whole contour.

    :param int nedges: number of line segments of resulting contour.
         Should be equal or greater than the number of sections defined by
         **pnts**.

    :returns: contour identifier
    """
    icheck(0, List(Point2D(), minlen=3))
    icheck(1, ListOr1(ZType(), llen=len(pnts) - 1))
    icheck(2, UInt(minv=len(pnts) - 1))

    b = bnds if isinstance(bnds, list) else [bnds]
    c = com.contcom.CreateSpline({"points": pnts,
                                  "bnds": b,
                                  "nedges": nedges})
    flow.exec_command(c)
    return c.added_contours2()[0]


@hmscriptfun
def extract_subcontours(source, plist, project_to="vertex"):
    """ Extracts singly connected subcontours from given contour.

        :param source: source contour identifier

        :param list-of-list-of-float plist: consecutive list of
           subcontours end points

        :param str project_to: defines projection rule for **plist** entries

           * ``"line"`` projects to closest point on the source contour
           * ``"vertex"`` projects to closest contour vertex
           * ``"corner"`` projects to closest contour corner vertex

        :returns: list of new contours identifiers

        Length of **plist** should be equal or greater than two.
        First and last points in **plist** define first and last points
        of **source** segment to extract.
        All internal points define internal division
        of this segment. Hence number of resulting subcontours will equal
        number of points in **plist** minus one.

        For closed **source** contour first and last **plist** points
        could coincide. In that case the sum of resulting subcontours
        will be equal to **source**.
    """
    icheck(0, ACont2D())
    icheck(1, List(Point2D(), minlen=2))
    icheck(2, OneOf("line", "vertex", "corner"))

    c = com.contcom.ExtractSubcontours({
        'src': source, 'plist': plist, 'project_to': project_to})
    flow.exec_command(c)
    return c.added_contours2()


@hmscriptfun
def add_rect_contour(p0, p1, bnd=0):
    """Adds four point closed rectangular contour.

    :param list-of-floats p0:

    :param list-of-floats p1: bottom left and top right coordinates of
       the contour

    :param bnd: single or list of 4 boundary
       identifiers (bottom, right, top, left) for contour segments.
       With the default value no boundary types will be set.

    :return: Contour identifier
    """
    icheck(0, Point2D())
    icheck(1, Point2D(grthan=p0))
    icheck(2, ListOr1(ZType(), llen=4))

    if isinstance(bnd, list):
        b = bnd[0:4]
    else:
        b = [bnd, bnd, bnd, bnd]

    c = com.contcom.AddRectCont({"p0": p0, "p1": p1, "bnds": b})
    flow.exec_command(c)
    return c.added_contours2()[0]


@hmscriptfun
def add_circ_contour(p0, rad, n_arc, bnd=0):
    """Adds circle contour from given center and radius.

    :param list-of-floats p0: circle center in [x, y] format

    :param float rad: circle radius

    :param int n_arc: partition of circle arc

    :param bnd: boundary identifier for contour.

    :return: Contour identifier
    """
    icheck(0, Point2D())
    icheck(1, Float(grthan=0.0))
    icheck(2, UInt(minv=3))
    icheck(3, ZType())

    c = com.contcom.AddCircCont({"p0": p0, "rad": rad,
                                 "na": n_arc, "bnd": bnd})
    flow.exec_command(c)
    return c.added_contours2()[0]


@hmscriptfun
def add_circ_contour2(p0, p1, p2, n_arc, bnd=0):
    """Adds circle contour from given arc points.

    :param list-of-floats p0:

    :param list-of-floats p1:

    :param list-of-floats p2: circle arc points as [x, y] format

    :param int n_arc: partition of circle arc

    :param bnd: boundary identifier for contour.
       With the default value no boundary types will be set.

    :return:  Contour identifier

    """
    icheck(0, Point2D())
    icheck(1, Point2D(noteq=[p0]))
    icheck(2, Point2D(noteq=[p0, p1]))
    icheck(3, UInt(minv=3))
    icheck(4, ZType())

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

    return add_circ_contour([cx + p0[0], cy + p0[1]], rad, n_arc, bnd)


@hmscriptfun
def add_circ_contour3(p0, p1, curv, n_arc, bnd=0):
    """Adds circle contour from given arc points and curvature.

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
    from hybmeshpack.basic.geom import angle_3pnt
    icheck(0, Point2D())
    icheck(1, Point2D(noteq=[p0]))
    icheck(2, Float(grthan=0.0))
    icheck(3, UInt(minv=3))
    icheck(4, ZType())

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
    a1 = angle_3pnt((0.0, 0.0), (cx1, cx2), (xa, ya))
    if (a1 < math.pi):
        cx, cy = cx1, cy1
    else:
        cx, cy = cx2, cy2
    return add_circ_contour([cx + p0[0], cy + p0[1]], r, n_arc, bnd)
