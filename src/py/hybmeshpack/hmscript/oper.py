from hybmeshpack import com
from hybmeshpack.hmscript import flow, data
from hybmeshpack.gdata.contour2 import Contour2
from hybmeshpack.basic.geom import Point2
from . import ExecError
import copy


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
       ValueError, hmscript.ExecError

    .. note::

       All contours from `cont` list are excluded consecutively.
       If you want to exclude multiply connected domain area you should first
       assemble multiply connected domain from the list of singly
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


def unite_grids(base_grid, over_grids, empty_holes=False, fix_bnd=False,
                zero_angle_approx=0):
    """Makes grids superpositions.

    :param base_grid: basic grid identifier

    :param list-of-tuples over_grids: sequence of grids for superposition as
      ``[(grid_id, buffer), () ...]`` where
      ``grid_id`` is an superposed grid identifier,
      ``buffer`` - size of the buffer for current imposition

    :param bool empty_holes: keep all empty zones
      (in case of multiple connectivity)
      of imposed grids in the resulting grid.

    :param bool fix_bnd: whether to fix all boundary nodes

    :param positive-degree zero_angle_approx:
      defines deviation from the straight angle which is considered
      insignificant. Grid boundary vertices which provide insignificant
      contour turns could be moved in order to obtain better result.
      Makes sense only if ``fix_bnd = False``.

    :return: identifier of the newly created grid

    :raises: hmscript.ExecError

    See detailed options description in :ref:`gridimp`.

    Example:
       .. literalinclude:: ../../testing/py/fromdoc/ex_unite_grids.py
           :start-after: START OF EXAMPLE
           :end-before: END OF EXAMPLE

    """
    args = {"base": base_grid, "empty_holes": empty_holes,
            "angle0": zero_angle_approx,
            "fix_bnd": fix_bnd, "plus": []}
    for ig in over_grids:
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
    :ivar str bnd_stepping: algorithm for stepping along the contour

      * ``'no'``: no artificial stepping. Only contour vertices will
        be used as grid nodes
      * ``'const'``: use only artificial stepping and ignore contour
        vertices
      * ``'keep_shape'``: use stepping and keep significant
        contour vertices
      * ``'keep_all'``: use stepping and keep all contour vertices
      * ``'incremental'``: increase boundary step linearly from ``start_point``
        to ``end_point``. Valid only for open contours
        ``start_point`` != ``end_point``. It acts like ``const`` stepping:
        no initial boundary vertices will be saved.

      .. note::

        While using ``const`` boundary step resulting grid may contain
        boundary edges with intermediate nodes. They are placed there
        intentionally to make further grid union operations easier.
        To get rid of such nodes use :func:`heal_grid` with
        ``simplify_boundary`` option at the very end of grid creation
        (see :ref:`example3` for example of elimination of those nodes).


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
    :ivar str project_to:
      option which defines **start_point**, **end_point** projection algorithm:

      * ``"line"`` - projects point to source contour line,
      * ``"vertex"`` - projects to closest source contour vertex,
      * ``"corner"`` - projects point to closest corner vertex.
    """

    def __init__(self, contour_id=None,
                 partition=[0],
                 direction="left",
                 bnd_step=0.1,
                 bnd_stepping="keep_shape",
                 range_angles=[40, 125, 235, 275],
                 force_conformal=False,
                 start_point=None,
                 end_point=None,
                 project_to="line"):
        """ Constructor with default attributes values given
        """
        self.contour_id = contour_id
        self.partition = copy.deepcopy(partition)
        self.direction = direction
        self.bnd_stepping = bnd_stepping
        self.bnd_step = bnd_step
        self.range_angles = copy.deepcopy(range_angles)
        self.force_conformal = force_conformal
        self.start_point = start_point
        self.end_point = end_point
        self.project_to = project_to

    def uniform_partition(self, fullh, n):
        """Sets uniform boundary grid vertical `partition` with
        the full depth of ``fullh`` and ``n`` total layers
        """
        h = fullh / n
        self.partition = [i * h for i in range(n + 1)]

    def incremental_partition(self, h0, coef, n):
        """Sets incremental boundary grid vertical `partition`
        starting from ``h0``, with each next step is ``coef`` times
        bigger then the previous. Total of ``n`` layers will be built

        """
        self.partition = [0, h0]
        for i in range(n - 1):
            self.partition.append(
                self.partition[-1] +
                coef * (self.partition[-1] - self.partition[-2]))


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
            if op.project_to == "line":
                pass
            elif op.project_to == "vertex" or op.project_to == "corner":
                c1 = data.get_any_contour(op.contour_id)
                if op.project_to == "corner":
                    c2 = Contour2.create_from_abstract(c1).simplify(0)
                    if c2 is not None:
                        c1 = c2
                i1 = c1.closest_point_index(d['start'])
                d['start'] = copy.deepcopy(c1.points[i1])
                i2 = c1.closest_point_index(d['end'])
                d['end'] = copy.deepcopy(c1.points[i2])
            else:
                raise ValueError("Unknown `project_to` = %s" % op.project_to)

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


def map_grid(base_grid, target_contour, base_points, target_points,
             snap="no", project_to="line", btypes="from_grid",
             algo="inverse_laplace", return_invalid=False):
    """Performs mapping of base grid on another contour.
    See detailed options description in :ref:`gridmappings`.

    :param base_grid: grid identifier.

    :param target_contour: contour identifier.

    :param base_points: collection of points in ``[[x0, y0], [x1, y1], ...]``
      format which lie on the **base_grid** contour (if a point doesn't lie on
      contour it would be projected to it).

    :param target_points: collection of points in ``[[x0, y0], [x1, y1], ...]``
      format which lie on the **target_contour**.
      The i-th point of **target_points** will be mapped into i-th point of
      **contour_points**. 

    :param str snap:
      an option which defines post processing algorithm of snapping
      newly created grid to *target_contour*:

      * ``"no"`` - no snapping
      * ``"add_vertices"`` - snap by adding new vertices if that will not
        ruin grid topology
      * ``"shift_vertices"`` - shift non-corner boundary nodes to corner
        locations if possible

    :param str project_to:
      option which defines **target_points** and **base_points**
      projection algorithm:

      * ``"line"`` - projects point to source contour line,
      * ``"vertex"`` - projects to closest source contour vertex,
      * ``"corner"`` - projects point to closest corner vertex.

    :param str btypes:
       defines from what source boundary features for newly created grid
       would be taken: ``"from_grid"`` or ``"from_contour"``.

    :param str algo:
       defines algorithm of mapping:

       * ``"direct_laplace"`` solves Laplace problem in base domain,
       * ``"inverse_laplace"`` solves Laplace problem in target domain.

    :param bool return_invalid: if this flag is on
       then the procedure will return a grid even if it is not valid
       (has self-intersections). Such grids could be exported to
       simple formats (like vtk or tecplot) in order to detect
       bad regions and give user a hint of how to adopt
       input data to gain an acceptable result.

       .. warning:: Never use invalid grids for further operations.

    :raises: ValueError, hmscript.ExecError

    :returns: identifier of newly created grid

    """
    n = max(len(base_points), len(target_points))
    bpoints = []
    tpoints = []
    for i in range(n):
        bpoints.append(Point2(*base_points[i]))
        tpoints.append(Point2(*target_points[i]))
    if project_to == "line":
        pass
    elif project_to == "vertex" or project_to == "corner":
        cbase = data.get_any_contour(base_grid)
        ctar = data.get_any_contour(target_contour)
        if project_to == "corner":
            cbase2 = Contour2.create_from_abstract(cbase).simplify(0)
            ctar2 = Contour2.create_from_abstract(ctar).simplify(0)
            if cbase2 is not None:
                cbase = cbase2
            if ctar2 is not None:
                ctar = ctar2
        for p in bpoints:
            i = cbase.closest_point_index(p)
            p.x, p.y = cbase.points[i].x, cbase.points[i].y
        for p in tpoints:
            i = ctar.closest_point_index(p)
            p.x, p.y = ctar.points[i].x, ctar.points[i].y
    else:
        raise ValueError("Unknown `project_to` = %s" % project_to)
    if algo not in ['inverse_laplace', 'direct_laplace']:
        raise ValueError("Unknown mapping algorithm")
    c = com.gridcom.MapGrid({"base": base_grid,
                             "target": target_contour,
                             "base_points": bpoints,
                             "target_points": tpoints,
                             "snap": snap,
                             "algo": algo,
                             "btypes": btypes,
                             "return_invalid": return_invalid})
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception as e:
        raise ExecError("map_grid. " + str(e))
