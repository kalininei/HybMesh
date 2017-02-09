"contours and 2D domains operations and modifications"
import copy
from hybmeshpack import com
from hybmeshpack.hmscript import flow, ExecError


def simplify_contour(cont, simplify=True, angle=0., separate=False):
    """ Separates and simplify user contour

    :param cont: source contour or grid identifier

    :param bool simplify: do simplification, i.e. make all
       segments non-collinear
       Collinear segments will not be splitted if they have different
       boundary types.

    :param degree angle: minimum allowed angle between segments
       after simplification (deg, >=0)

    :param bool separate: assemble list of singly connected contours from
       multiply connected source contour

    :returns: list of created contours ids

    """
    c = com.contcom.SimplifyContours({"names": [cont],
                                      "simplify": simplify,
                                      "angle": angle,
                                      "separate": separate})
    try:
        flow.exec_command(c)
        return c.added_contours2()
    except:
        raise ExecError("simplify_contour")


def decompose_contour(cont):
    """ Returns set of simple singly connected contours built from
    input contour

    :param cont: contour identifier

    :returns: list of new contour identifiers

    :raises: ValueError, ExecError

    All resulting contour vertices will have no more than two adjacent edges.
    If input contour has vertices with more than two connections
    (or self intersections) it will be splitted at these points.

    Tries to assemble as many closed contours as possible.
    """
    c = com.contcom.DecomposeContour({"source": cont})
    try:
        flow.exec_command(c)
        return c.added_contours2()
    except:
        raise ExecError("decompose_contour")


def unite_contours(conts):
    """ Unites contours to single multiply connected contour

    :param conts: list of contours identifiers to unite

    :return: contour identifier

    Equal nodes of input contours will be merged.
    """
    c = com.contcom.UniteContours({"sources": conts})
    try:
        flow.exec_command(c)
        return c.added_contours2()[0]
    except Exception as e:
        raise ExecError("unite_contours. " + str(e))


def set_boundary_type(cid, btps=None, bfun=None):
    """ Mark user or grid contour segments with boundary types.

    :param cid: contour or grid identifier

    :param btps: list of boundary types identifiers for each segment
       or single identifier for the whole contour

    :param  bfun: function ``(x0, y0, x1, y1, bt) -> btype``
       which returns boundary type taking segment endpoints
       and old boundary type as arguments.

    Example:

      .. literalinclude:: ../../testing/py/fromdoc/ex_setbtype.py
          :start-after: START OF EXAMPLE
          :end-before: END OF EXAMPLE

    Only one of **btps**, **bfun** arguments should be defined.
    """
    try:
        cont = flow.receiver.get_any_contour(cid)
        if bfun is not None:
            bt = cont.raw_data('btypes')
            pt = cont.raw_data('vertices')
            ev = cont.raw_data('edge-vert')
            newbc = []
            for ibt, iev in zip(bt, ev):
                p0, p1 = pt[iev[0]], pt[iev[1]]
                newbc.append(bfun(p0[0], p0[1], p1[0], p1[1], ibt))
            for i in range(len(newbc)):
                if newbc[i] is None:
                    newbc[i] = bt[i]
        else:
            newbc = btps
        args = {}
        if not isinstance(newbc, list):
            args[newbc] = range(cont.n_edges())
        else:
            for i, b in enumerate(newbc):
                if b not in args:
                    args[b] = []
                args[b].append(i)

        c = com.contcom.SetBTypeToContour({"name": cid, "btypes": args})
        flow.exec_command(c)
    except:
        raise ExecError("set_boundary_type")


def clip_domain(dom1, dom2, operation, simplify=True):
    """ Executes domain clipping procedure

    :param dom1:

    :param dom2: contour identifiers

    :param  str operatrion: operation code

       * ``"union"``
       * ``"difference"``
       * ``"intersection"``
       * ``"xor"``

    :param bool simplify: whether to keep all source points (False) or
       return simplified contour

    :returns: created contour identifier or None if resulting domain is empty
    """
    c = com.contcom.ClipDomain({"c1": dom1, "c2": dom2, "oper": operation,
                                "simplify": simplify})
    try:
        flow.exec_command(c)
        if len(c.added_contours2()) > 0:
            return c.added_contours2()[0]
        else:
            return None
    except:
        raise ExecError("clip_domain")


def partition_contour(cont, algo, step=1., angle0=30., keep_bnd=False,
                      nedges=None, crosses=[], keep_pts=[],
                      start=None, end=None):
    """ Makes connected contour partition

    :param cont: Contour or grid identifier

    :param str algo: Partition algorithm:

       * ``'const'``: partition with defined constant step
       * ``'ref_points'``: partition with step function given by a
         set of values refered to basic points
       * ``'ref_weights'``: partition with step function given by a
         set of values refered to local contour [0, 1] coordinate
       * ``'ref_lengths'``: partition with step function given by a
         set of values refered to local contour length coordinate

    :param step: For ``algo='const'`` a float number defining
       partition step;

       For ``algo='ref_points'`` - list of step values and point coordinates
       given as
       ``[ step_0, [x_0, y_0], step_1, [x_1, y_1], ....]``.

       For ``algo='ref_weights'`` - list of step values and point normalized
       1d coordinates given as
       ``[ step_0, w_0, step_1, w_1, ....]``.

       For ``algo='ref_lengths'`` - list of step values and point 1d
       coordinates given as
       ``[ step_0, s_0, step_1, s_1, ....]``. Negative ``s_i`` shows
       length coordinate started from end point in reversed direction.


    :param float angle0: existing contour vertices which provide
       turns outside of ``[180 - angle0, 180 + angle0]`` degrees range
       will be preserved regardless of other options

    :param bool keep_bnd: if that is True than vertices which have different
       boundary features on their right and left sides will be preserved

    :param int nedges: if this parameter is not None then it provides
       exact number of edges in the resulting contour. To satisfy this
       condition **step** value will be multiplied by an appropriate factor.
       If it can not be satisfied (due to other restrictions)
       then an exception will be raised.

    :param crosses: represents set of contour which cross points with
       target contour will be present in resulting contour.

    :param keep_pts: list of points as ``[[x1, y1], [x2, y2], ...]`` list
       which should present in output partition.

    :param start:

    :param end: start and end points which define processing segment

    :returns: new contour identifier

    :raises: hmscript.ExecError, ValueError

    Points set defined by user for ``algo='ref_points'`` algorithm
    will not present in resulting contour (as well as points defined
    implicitly by other ``'ref_'`` algorithms). It just shows locations
    where step size of given length should be applied. If any point
    of this set is not located on the input contour then it will be
    projected to it.

    For constant stepping any contour including multiply connected ones
    could be passed. For ``ref_`` stepping only singly connected
    contours (open or closed) are allowed.

    If **start** and **end** points are defined then only segment between
    these points will be parted. Note that all closed contours are treated
    in counterclockwise direction. Given start/end points will be
    projected to closest contour vertex.

    ``ref_weights``, ``ref_lengths`` partition methods
    require definition of **start** point
    (To process whole contour ``end`` point could be omitted).

    Example:

      .. literalinclude:: ../../testing/py/fromdoc/ex_partcontour.py
          :start-after: vvvvvvvvvvvvvvvvvvvvvvvv
          :end-before: ^^^^^^^^^^^^^^^^^^^^^^^^

    See also: :ref:`simplecontmeshing`
    """
    # prepare arguments for command
    if algo == "const":
        plain_step = [step]
    elif algo == "ref_points":
        plain_step = []
        for i in range(len(step) / 2):
            plain_step.append(step[2 * i])
            plain_step.extend(step[2 * i + 1])
    elif algo in ["ref_weights", "ref_lengths"]:
        plain_step = copy.deepcopy(step)
    else:
        raise ValueError
    sp, ep = None, None
    if start is not None:
        sp = [start[0], start[1]]
    if end is not None:
        ep = [end[0], end[1]]
    kp = []
    if keep_pts is not None:
        for p in keep_pts:
            kp.append([p[0], p[1]])

    args = {"algo": algo,
            "step": plain_step,
            "angle0": angle0,
            "keepbnd": keep_bnd,
            "base": cont,
            "nedges": nedges,
            "crosses": crosses,
            "start": sp,
            "end": ep,
            "keep_pts": kp}
    # call
    c = com.contcom.PartitionContour(args)
    try:
        flow.exec_command(c)
        return c.added_contours2()[0]
    except Exception:
        raise ExecError('partition_contour')


def matched_partition(cont, step, influence, ref_conts=[], ref_pts=[],
                      angle0=30., power=3.):
    """ Makes a contour partition with respect to other
        contours partitions and given reference points

        :param cont: target contour.

        :param float step: default contour step size.

        :param float influence: influence radius of size conditions.

        :param ref_conts: list of contours which segmentation will
          be treated as target segmentation conditions.

        :param ref_pts: reference points given as
           ``[step0, [x0, y0], step1, [x1, y1], ...]`` list.

        :param float angle0: existing contour vertices which provide
           turns outside of ``[180 - angle0, 180 + angle0]`` degrees range
           will be preserved regardless of other options

        :param positive-float power: shows power of weight
           calculation function. As this parameter increases
           size transitions along contour become less smooth and more
           sensible to size conditions.

        See :ref:`matchedcontmeshing` for details.
    """
    if not isinstance(ref_conts, list):
        ref_conts = [ref_conts]

    rp = []
    for i in range(len(ref_pts) / 2):
        rp.append(ref_pts[2 * i])
        rp.extend(ref_pts[2 * i + 1])

    args = {"base": cont,
            "cconts": ref_conts,
            "cpts": rp,
            "step": step,
            "infdist": influence,
            "angle0": angle0,
            "power": power
            }
    c = com.contcom.MatchedPartition(args)
    try:
        flow.exec_command(c)
        return c.added_contours2()[0]
    except Exception:
        raise ExecError('matched_partition')


def partition_segment(start, end, hstart, hend, hinternal=[]):
    """ Makes a partition of numeric segment by given
    recommended step sizes at different locations.

    :param float start:

    :param float end: start and end points of numeric segment.

    :param float hstart:

    :param float hend: recommended step at the beginning and at the
       end of numeric segments

    :param list-of-floats hinternal: possible internal recommended steps
       given as [segment value0, step0, segment value1, step1, ...]

    :returns: increasing float list representing
       partitioned segment [start, ..(internal steps).., end]

    :raises: ValueError, ExecError

    This is a helper function which runs iterative procedure
    to adopt partition of the segment to user defined step size
    recommendations. It could be used for example to
    explicitly define partition in :func:`add_unf_rect_grid`,
    :func:`add_unf_circ_grid` with given `custom\_\*` options or
    vertical partition of boundary grid in :func:`build_boundary_grid`.

    Example:
      Shows division of [1, 2] segment from approximately 0.1
      step size at the beginning to approximately 0.5 step size at
      the end

      .. code-block:: python

         print hm.partition_segment(1, 2, 0.1, 0.5)
         >>> [1.0, 1.1238348, 1.309012, 1.585923, 2.0]

      Shows division of [1, 2] segment with approximate 0.1
      step sizes at its end points and 0.4 step size at the center.

      .. code-block:: python

         print hm.partition_segment(1, 2, 0.1, 0.1, [1.5, 0.4])
         >>> [1.0, 1.1082238, 1.3278488, 1.672151, 1.891776, 2.0]
    """
    from hybmeshpack.hmcore import c2 as c2core
    try:
        return c2core.segment_partition(start, end, hstart, hend, hinternal)
    except:
        raise ExecError("partition_segment")


def connect_subcontours(sources, fix=[], close="no", shiftnext=True):
    """ Connects sequence of open contours into a single contour
    even if neighboring contours have no equal end points

    :param sources: list of open contour identifiers

    :param list-of-int fix: indicies of **sources** contours
        which could not be shifted or stretched.

    :param str close: last connection algorithm:
        ``no``, ``yes`` or ``force``

    :param bool shiftnext: if True then each next contour will be
        shifted to the end point of previous one, otherwise
        both contours will be stretched to match average end point.

    To connect given contours this procedure implements stretching
    and shifting of those ones not listed in **fix** list.
    If two adjacent source contours are marked as fixed but have no
    common end points an exception will be raised.

    If **close** is ``yes`` then last contour of **sources** will
    be connected with the first one with algorithm depending on
    **fix** and **shiftnext** options.

    If **close** is ``no``
    then ending contours will be left as they are. In that case resulting
    contour will be open until first and last points are exactly equal.

    **close** = ``force`` algorithm works like **close** = ``no`` but
    creates a section which explicitly connects first and last contours
    by a straight line.
    """
    args = {}
    args['src'] = sources
    args['fix'] = copy.deepcopy(fix)
    args['close'] = close
    args['shiftnext'] = shiftnext
    c = com.contcom.ConnectSubcontours(args)
    try:
        flow.exec_command(c)
        return c.added_contours2()[0]
    except Exception:
        raise ExecError("connect_subcontours")
