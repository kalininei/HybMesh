from HybMeshPyPack import com
from HybMeshPyPack.basic.geom import Point2
from HybMeshPyPack.hmscript import flow, data


def add_rect_cont(p0, p1, bnd=None):
    """Adds four point rectangular contour

    Args:
       p0, p1 (list-of-floats): bottom left and top right coordinates of
       the contour

    Kwargs:
       bnd (boundary identifier): single or list of 4 boundary
       identifiers (bottom, left, top, right) for contour segments.
       With the default value no boundary types will be set.

    Returns:
       Contour identifier

    Example:
       >>> hmscript.AddRectCont([0, 0], [1,1], [b1, b2, b1, b2])

    """
    if isinstance(bnd, list):
        b = bnd[0:4]
    elif bnd is not None:
        b = [bnd, bnd, bnd, bnd]
    else:
        b = [0, 0, 0, 0]
    c = com.contcom.AddRectCont({"p0": Point2(*p0),
                                 "p1": Point2(*p1),
                                 "bnds": b})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def grid_bnd_to_cont(g1, simplify=True):
    """ Extracts grid boundary to user contour

    Args:
       g1: grid identifier

    Kwargs:
       simplify (bool): if true deletes all non-significant points. Otherwise
       resulting contour will contain all boundary grid edges.

    Returns:
       contour identifier
    """
    c = com.contcom.GridBndToContour({"grid_name": g1,
                                      "simplify": simplify})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def simplify_contour(cont, simplify=False, angle=0, separate=False):
    """ Separates and simplify user contour

    Args:
       cont: source contour or grid identifier

    Kwargs:
       simplify (bool): do simplification, i.e. make all segments non-collinear
       Collinear segments will not be splitted if they have different
       boundary types.

       angle (float): minimum acute angle between segments
       after simplification (deg, >=0)

       separate (bool): assemble list of singly connected contours from
       multiply connected source contour

    Returns:
        list of created contours ids

    """
    c = com.contcom.SimplifyContours({"names": [cont],
                                      "simplify": simplify,
                                      "angle": angle,
                                      "separate": separate})
    flow.exec_command(c)
    return c._get_added_names()[1]


def unite_contours(conts):
    """ Unites contours to single multiply-connected contour

    Args:
       conts: list of contours identifiers to unite

    Returns:
       contour identifier

    .. note::

       Input contours should not overlap each other.
       This procedure doesn't make checks for it.

    """
    c = com.contcom.UniteContours({"sources": conts})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def add_boundary_type(index, name="boundary1"):
    """ Creates boundary type

    Args:
       index (int) - index of boundary (>0)

    Kwargs:
       name (str) - user defined name of the boundary

    Returns:
       boundary identifier

    If boundary with `index` already exists it will be overwriten.
    Name of the boundary should be unique. If it already exists it will
    be changed.

    """
    c = com.contcom.EditBoundaryType({"index": index, "name": name})
    flow.exec_command(c)
    return index


def set_btype_to_contour(cont, btps=None, bfun=None):
    """ Mark user or grid contour segments with btypes.

    Args:
       cont - contour or grid identifier

    Kwargs:
       btps - list of boundary types identifiers for each segment
       or single identifier for the whole contour

       bfun - function (x0, y0, x1, y1, bt) -> btype
       which returns boundary type from segment endpoints and old boundary type

    Example:

      .. code-block:: python

         # create rectangular contour
         cont = hmscript.AddRectContour([0, 0], [1, 1])
         # create boundary types
         b1 = hmscript.AddBoundaryType(1, "vertical")
         b2 = hmscript.AddBoundaryType(2, "horizontal")
         # mark contour using `btps`
         hmscript.SetBTypeToContour(cont, btps=[b2, b1, b2, b1])
         # same command using `bfun`
         def bfun(x0, y0, x1, y1, bold):
            global b1, b2
            return b1 if x0 == x1 else b2
         hmscript.SetBTypeToContour(cont, bfun=bf)

    Only one of `btps`, `bfun` arguments should be defined.
    """
    cdata = data.get_any_contour(cont)
    args = {"name": cont}
    if btps is not None:
        if not isinstance(btps, list):
            args[btps] = range(cdata.n_edges())
        else:
            for i, b in enumerate(btps):
                if b not in args:
                    args[b] = []
                args[b].append(i)
    elif bfun is not None:
        ep = cdata.edges_points()
        for i, iep in enumerate(ep):
            p0 = cdata.points[iep[0]]
            p1 = cdata.points[iep[1]]
            r = bfun(p0.x, p0.y, p1.x, p1.y, cdata.edge_bnd(i))
            if r is None:
                continue
            if r not in args:
                args[r] = []
            args[r].append(i)

    c = com.contcom.SetBTypeToContour({"btypes": [args]})
    flow.exec_command(c)


def create_contour(pnts, bnds=0):
    """ Create singly connected contour from sequence of points

    Args:
       pnts(list-of-list-of-floats): sequence of points

    Kwargs:
       bnds(single or list-of boundary identifiers): boundary type for
       each contour segment or single identifier for the whole contour.

    Example:
       >>> hmscript.CreateContour([[0, 0], [1, 0], [1, 1], [0, 0]])

    If coordinates of first and last points are equal
    then contour is considered closed.

    """
    pt = []
    for p in pnts:
        pt.append(Point2(*p))
    c = com.contcom.CreateContour({"points": pt,
                                   "bnds": bnds})
    flow.exec_command(c)
    return c._get_added_names()[1][0]
