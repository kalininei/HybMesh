from hybmeshpack import com
from hybmeshpack.basic.geom import Point2
from hybmeshpack.hmscript import flow, data
from hybmeshpack.hmscript import ExecError


def grid_bnd_to_contour(g1, simplify=True):
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


def simplify_contour(cont, simplify=True, angle=0, separate=False):
    """ Separates and simplify user contour

    Args:
       cont: source contour or grid identifier

    Kwargs:
       simplify (bool): do simplification, i.e. make all segments non-collinear
       Collinear segments will not be splitted if they have different
       boundary types.

       angle (float): minimum allowed angle between segments
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
    """ Unites contours to single multiply connected contour

    Args:
       conts: list of contours identifiers to unite

    Returns:
       contour identifier

    .. note::

       Input contours should not cross each other.
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
       integer boundary identifier

    If boundary with ``index`` already exists it will be overwriten.
    Name of the boundary should be unique, if it already exists it will
    be changed automatically.

    """
    c = com.contcom.EditBoundaryType({"index": index, "name": name})
    flow.exec_command(c)
    return index


def set_boundary_type(cont, btps=None, bfun=None):
    """ Mark user or grid contour segments with boundary types.

    Args:
       cont - contour or grid identifier

    Kwargs:
       btps - list of boundary types identifiers for each segment
       or single identifier for the whole contour

       bfun - function (x0, y0, x1, y1, bt) -> btype
       which returns boundary type from segment endpoints
       and old boundary type

    Example:

      .. literalinclude:: ../../testing/py/fromdoc/ex_setbtype.py
          :start-after: START OF EXAMPLE
          :end-before: END OF EXAMPLE

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
       pnts(list-of-list-of-floats): sequence of points.
       If coordinates of first and last points are equal
       then contour is considered closed.

    Kwargs:
       bnds(single or list-of boundary identifiers): boundary type for
       each contour segment or single identifier for the whole contour.

    Example:
       >>> hmscript.create_contour([[0, 0], [1, 0], [1, 1], [0, 0]],
                                  [b1, b2, b3])
    """
    pt = []
    for p in pnts:
        pt.append(Point2(*p))
    c = com.contcom.CreateContour({"points": pt,
                                   "bnds": bnds})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def clip_domain(dom1, dom2, operation, simplify=True):
    """ Executes domain clipping procedure

    Args:
       dom1, dom2: contour identifiers

       operatrion (str): operation code
          * ``"union"``
          * ``"difference"``
          * ``"intersection"``
          * ``"xor"``

    Kwargs:
       simplify (bool): whether to keep all source points (False) or
       return simplified contour

    Returns:
       created contour identifier or None if resulting domain is empty
    """
    if operation not in ['union', 'difference', 'intersection', 'xor']:
        raise ValueError("unknows operation: %s" % str(operation))
    c = com.contcom.ClipDomain({"c1": dom1, "c2": dom2, "oper": operation,
                                "simplify": simplify})
    flow.exec_command(c)
    if len(c._get_added_names()[1]) > 0:
        return c._get_added_names()[1][0]
    else:
        return None


def partition_contour(cont, algo, step, angle0=30, keep_bnd=False):
    """ Makes singly connected contour partition

    :param cont: Contour or grid identifier

    :param algo: Partition algorithm:
       *``'const'``
       *``'ref_points'``

    :param step: For ``algo="const"`` a float number defining
       partition step;

       For ``algo='ref_points'`` - list of step value and point coordinates:
       ``[ step_0, [x_0, y_0], step_1, [x_1, y_1], ....]``

    :param zero_angle_approx: existing contour vertices which provide
       turns in ``[180 - angle0, 180 + angle0]`` degrees range are considered
       insignificant and will not be preserved

    :param keep_bnd: if that is True than vertices which have different
       boundary features on their right and left sides will be preserved

    :returns: new contour identifier

    :raises: hmscript.ExecError, ValueError
    """
    import numbers
    # checks
    if algo == "const":
        if not isinstance(step, numbers.Real):
            raise ValueError("invalid step for const contour partition")
    elif algo == "ref_points":
        errs = "invalid step for ref_points partition"
        if not isinstance(step, list) or len(step) % 2 != 0:
            raise ValueError(errs)
        for i in range(len(step) / 2):
            if not isinstance(step[2 * i], numbers.Real):
                raise ValueError(errs)
            if not isinstance(step[2 * i + 1], list) or\
                    len(step[2 * i + 1]) != 2 or\
                    not isinstance(step[2 * i + 1][0], numbers.Real) or\
                    not isinstance(step[2 * i + 1][1], numbers.Real):
                raise ValueError(errs)
    else:
        raise ValueError("unknown partition angorithm")
    # prepare arguments for command
    if algo == "const":
        plain_step = [step]
    elif algo == "ref_points":
        plain_step = []
        for i in range(len(step) / 2):
            plain_step.append(step[2 * i])
            plain_step.extend(step[2 * i + 1])

    args = {"algo": algo,
            "step": plain_step,
            "angle0": angle0,
            "keepbnd": keep_bnd,
            "base": cont}
    # call
    c = com.contcom.PartitionContour(args)
    try:
        flow.exec_command(c)
        return c._get_added_names()[1][0]
    except Exception:
        raise ExecError('partition contour')
