from hybmeshpack import com
from hybmeshpack.hmscript import flow
from hybmeshpack.basic.geom import Point2
from . import ExecError


def exclude_contours(grid, conts, what):
    """Builds a grid by excluding contour area from existing grid.

    Args:
       grid: source grid identifier

       conts: contour or list of contours/grids identifiers for exclusion.

       what (str): ``"inner"``/``"outer"``.
       Describes what part of ``conts`` domain exclude

    Returns:
       new grid identifier

    Raises:
       hmscript.ValueError, hmscript.ExecError

    .. note::

       All countours from `cont` list are excluded consecutively.
       If you want to exclude multiply connected domain area you should first
       assemle multiply connected domain from the list of singly
       connected ones using
       :func:`unite_contours` procedure.

    Example:
       .. literalinclude:: ../../testing/py/fromdoc/ex_exclude.py
           :start-after: START OF EXAMPLE
           :end-before: END OF EXAMPLE

    """
    if not isinstance(conts, list):
        conts = [conts]
    if what == "inner":
        exa = True
    elif what == "outer":
        exa = False
    else:
        raise ValueError("Invalid domain area: %s" % str(what))
    c = com.gridcom.ExcludeContours({"grid_name": grid,
                                     "cont_names": conts,
                                     "is_inner": exa})
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except:
        raise ExecError('exclude_contours')


def unite_grids(base_grid, imp_grids, empty_holes=False, fix_bnd=False):
    """Makes grids impositions

    Args:
      base_grid: basic grid identifier

      imp_grids (list-of-tuples): sequence of grids for imposition as
      ``[(grid_id, buffer), () ...]`` where
      ``grid_id`` is an imposed grid identifier,
      ``buffer`` - size of the buffer for current imposition

    Kwargs:
      empty_holes (bool): keep all empty zones (in case of multle connectivity)
      of imposed grids in the resulting grid.

      fix_bnd (bool): whether to fix all boundary nodes

    Returns:
       identifier of the newly created grid

    Raises:
       hmscript.ExecError

    Each next grid will be imposed on the result of previous imposition.

    Example:
       .. literalinclude:: ../../testing/py/fromdoc/ex_unite_grids.py
           :start-after: START OF EXAMPLE
           :end-before: END OF EXAMPLE

    """
    args = {"base": base_grid, "empty_holes": empty_holes,
            "fix_bnd": fix_bnd, "plus": []}
    for ig in imp_grids:
        args["plus"].append({"name": ig[0], "buf": ig[1], "den": 7})
    c = com.gridcom.UniteGrids(args)
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('unite_grids')


class BoundaryGridOptions(object):
    """Options for building boundary grid

    :ivar contour_id:
      identifier of input source contour
    :ivar list-of-floats partition:
      partition in perpendicular direction.
      List of ascending floats starting from zero which
      represents distances from contour to grid layers
    :ivar str direction:
      'left'/'right'. Boundary grid will be built to the left/right
      from the contour (with respect to positive tracing)
    :ivar float bnd_step:
      float size of artificial stepping along the contour (used only if
      ``bnd_stepping`` is not "no"). If ``bnd_stepping == 'incremental'`` then
      this should be a list of two floats: boundary partition nearby start
      and end of the contour
    :ivar str bnd_stepping:
      algorithm for stepping along the contour

      * ``'no'``: no artificial stepping. Only contour vertices will
        be used as grid nodes
      * ``'const'``: use only artificial stepping and ignore contour
        vertices
      * ``'keep_shape'``: use stepping and keep significant
        contour vertices
      * ``'keep_all'``: use stepping and keep all contour vertices
      * ``'incremental'``: increase boundary step lineary from ``start_point``
        to ``end_point``. Valid only for open contours
        ``start_point`` != ``end_point``. It acts like ``const`` stepping:
        no initial boundary vertices will be saved.
    :ivar list-of-floats range_angles:
      list of 4 angle values (deg) which define algorithms
      for contour bends treatment:

      * ``[0,     ra[0]]``: acute angle algorithm
      * ``[ra[0], ra[1]]``: right angle algorithm
      * ``[ra[1], ra[2]]``: straight angle algorithm
      * ``[ra[2], ra[3]]``: reentrant angle algorithm
      * ``[ra[3],   360]``: round algorithm
    :ivar bool force_conformal:
      use strictly conformal mappings
    :ivar list-of-floats start_point, end_point:
      points in [x, y] format which define
      the exact segment of the contour for building grid.
      If both are None hence whole contour (or all subcontours) will be used.
      if ``start_point`` == ``end_point`` then whole subcontour closest to
      this point will be used.
      If point is not located on the source contour then it will be
      projected to it.
    """

    def __init__(self, contour_id=None,
                 partition=[0],
                 direction="left",
                 bnd_step=0.1,
                 bnd_stepping="keep_shape",
                 range_angles=[40, 125, 235, 275],
                 force_conformal=False,
                 start_point=None,
                 end_point=None):
        """ Constructor with default attributes values given
        """
        self.contour_id = contour_id
        self.partition = partition
        self.direction = direction
        self.bnd_stepping = bnd_stepping
        self.bnd_step = bnd_step
        self.range_angles = range_angles
        self.force_conformal = force_conformal
        self.start_point = start_point
        self.end_point = end_point

    def uniform_partition(self, fullh, n):
        """Sets uniform boundary grid vertical `partition` with
        the full depth of ``fullh`` and ``n`` total layers
        """
        h = fullh / n
        self.partition = [i * h for i in range(n + 1)]

    def incremental_partition(self, h0, coef, n):
        """Sets incremental boundary grid vertical `partition`
        starting from ``h0``, whith each next step is ``coef`` times
        bigger then the previous. Total of ``n`` layers will be built

        """
        self.partition = [0, h0]
        for i in range(n - 1):
            self.partition.append(
                self.partition[-1] +
                coef * (self.partition[-1] - self.partition[-2]))

    def reentrant_all_square(self):
        """This guaranties that all reentrant angles will be treated
        with additional square zone and without rounding.
        At angles close to 360 degree this may lead to a vary bad cells"""
        self.range_angles[3] = max(self.range_angles[3], 359)

    def reentrant_all_round(self):
        """This guaranties that all reentrant angles will be treated
        with round algorithm.
        """
        self.range_angles[3] = self.range_angles[2]


def build_boundary_grid(opts):
    """Builds a boundary grid near contour

    Args:
       opts: single or list of
       :class:`BoundaryGridOptions` objects

    Returns:
       identifier of the newly created grid.

    Raises:
       hmscript.ExecError, ValueError

    .. note::
      If different options for different segments of the contour are required
      (e.g. different partitions) then multiple option instances with the same
      target contour should be passed.

    Example:
       this code
       creates boundary grid along a square with different settings
       for vertical and horizontal edges of the source square:

       .. literalinclude:: ../../testing/py/fromdoc/ex_bgrid.py
           :start-after: START OF EXAMPLE
           :end-before: END OF EXAMPLE


    """
    inp = []
    if not isinstance(opts, list):
        opts2 = [opts]
    else:
        opts2 = opts
    for op in opts2:
        d = {}
        d['source'] = op.contour_id
        d['direction'] = op.direction
        if op.direction == "left":
            d['direction'] = 1
        elif op.direction == "right":
            d['direction'] = -1
        else:
            raise ValueError("Invalid direction: " + str(op.direction))
        if op.bnd_stepping == "no":
            d['mesh_cont'] = 0
        elif op.bnd_stepping == "const":
            d['mesh_cont'] = 3
        elif op.bnd_stepping == "keep_shape":
            d['mesh_cont'] = 2
        elif op.bnd_stepping == "keep_all":
            d['mesh_cont'] = 1
        elif op.bnd_stepping == "incremental":
            d['mesh_cont'] = 4
        else:
            raise ValueError("Invalid bnd_stepping: " + str(op.bnd_stepping))
        if op.bnd_stepping != "incremental":
            if isinstance(op.bnd_step, list):
                raise ValueError("list values for bnd_step are available "
                                 "only for incremental stepping")
            else:
                d['mesh_cont_step'] = op.bnd_step
        else:
            d['mesh_cont_step'] = 0.0
            if not isinstance(op.bnd_step, list) or len(op.bnd_step) < 2:
                raise ValueError("Incremental stepping requires list[2] as "
                                 "bnd_step")
            d['step_start'] = op.bnd_step[0]
            d['step_end'] = op.bnd_step[1]

        d['algo_acute'] = op.range_angles[0]
        d['algo_right'] = op.range_angles[1]
        d['algo_straight'] = op.range_angles[2]
        d['algo_reentr'] = op.range_angles[3]
        if (op.start_point is not None and op.end_point is not None):
            d['start'] = Point2(*op.start_point)
            d['end'] = Point2(*op.end_point)
        d['force_conf'] = op.force_conformal
        # check partition
        for i in range(len(op.partition) - 1):
            if (op.partition[i] >= op.partition[i + 1]):
                op.partition = [0]
                break
        if len(op.partition) < 2:
            raise ValueError("Invalid partition")
        d['partition'] = op.partition
        inp.append(d)

    c = com.gridcom.BuildBoundaryGrid({"opt": inp})
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('build_boundary_grid')
