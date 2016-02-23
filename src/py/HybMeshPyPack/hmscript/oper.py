from HybMeshPyPack import com
from HybMeshPyPack.hmscript import flow
from HybMeshPyPack.basic.geom import Point2
from . import ExecError


def exclude_contours(grid, conts, exclude_outer=True):
    """Builds a grid by excluding contour area from existing grid.

    Args:
       grid: source grid identifier

       conts: contour or list of contours/grids identifiers for exclusion.

    Kwargs:
       exclude_outer - exclude inner or outer region of contours

    Returns:
       new grid identifier

    Raises:
       hmscript.ExecError

    All input contours should bound a domain. Open contours are not allowed.

    .. note::

       All countours from `cont` list are excluded consecutively.
       If you want to exclude multiply connected domain area you should first
       assemle multiply connected domain by list of singly connected ones
       using :func:`HybMeshPyPack.hmscript.UniteContours` procedure.

    Example:

       .. code-block:: python

          # source grid
          g1 = hmscript.AddUnfCircGrid([0, 0], 5, 16, 5)
          # outer domain border
          c1 = hmscript.AddRectContour([-4, -4], [4, 4])
          # inner domain border
          c2 = hmscript.CreateContour([[-2, -2], [2, -2], [0, 2], [-2, -2]])
          # multiply connected domain
          c3 = hmscript.UniteContours([c1, c2])
          # exclusion
          g2 = hmscript.ExcludeContours(g2, c3)
    """
    if not isinstance(conts, list):
        conts = [conts]
    c = com.gridcom.ExcludeContours({"grid_name": grid,
                                     "cont_names": conts,
                                     "is_inner": not exclude_outer})
    flow.exec_command(c)
    return c._get_added_names()[0][0]


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

    Each next grid will be imposed on the result of previous imposition

    Example:

       .. code-block:: python

          # lower level grid
          g1 = hmscript.AddUnfRectGrid([0, 0], [10, 10], 10, 10)
          # first imposition grid
          g2 = hmscript.AddUnfRectGrid([0, 0], [3, 3], 7, 7)
          # second imposition grid
          g3 = hmscript.AddUnfCircGrid([5, 5], 3, 20, 8)
          # impose grids
          impgrid = hmscript.UniteGrids(g1, [(g1, 2.0), (g2, 2.0)])

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
        raise ExecError('build_boundary_grid')


class BoundaryGridOption(object):
    """Options for building boundary grid

    :ivar contour_id:
      identifier of input contour or grid
    :ivar list-of-floats partition:
      partition in perpendicular direction.
      List of ascending floats starting with zero which
      represents distances from contour to grid layers
    :ivar int direction:
      'left'/'right'. Boundary grid will be built to the left/right
      from the contour (with respect to positive tracing)
    :ivar float bnd_step:
      float size of artificial stepping along the contour (used only if
      `bnd_stepping` is not "no")
    :ivar str bnd_stepping:
      algorithm for stepping along the contour

      * 'no': no artificial stepping. Only contour edges will
        be used as grid nodes
      * 'const': use artificial stepping and ignore contour
        nodes
      * 'keep_shape': use stepping and keep significant
        contour nodes
      * 'keep_all': use stepping and keep all contour nodes
    :ivar list-of-floats range_angles:
      list of 4 angle (deg) values which define algorithms
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
      If both are None hence whole contour will be used.
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


def build_boundary_grid(opts):
    """Builds a boundary grid near contour

    Args:
       opts: single or list of
       :class:`HybMeshPyPack.hmscript.BoundaryGridOption` objects

    Returns:
       identifier of the newly created grid.

    Raises:
       hmscript.ExecError, ValueError

    .. note::
      If different options for different segements of the contour are required
      (e.g. different partitions) then multiple option instances with the same
      target contour should be passed. For example this code
      creates boundary grid along a square with different settings
      for vertical and horizontal edges of the source square:

      .. code-block:: python

        # create source contour
        cont = hmscript.AddRectCont([0, 0], [1, 1])
        # basic options for horizontal segments
        opvert = hmscript.BoundaryGridOption(cont,
           partition=[0, 0.01, 0.02, 0.03],
           bnd_stepping='keep_shape',
           bnd_step=0.05)
        # basic options for vertical segments
        ophoriz = hmscript.BoundaryGridOption(cont,
           partition=[0, 0.01, 0.02, 0.03, 0.05, 0.1],
           bnd_stepping='keep_shape',
           bnd_step=0.03)
        # option for bottom segment
        op1 = copy.deepcopy(ophoriz)
        op1.start_point = [0, 0]
        op1.end_point = [1, 0]
        # option for left segment
        op2 = copy.deepcopy(opvert)
        op2.start_point = [1, 0]
        op2.end_point = [1, 1]
        # option for top segment
        op3 = copy.deepcopy(ophoriz)
        op3.start_point = [1, 1]
        op3.end_point = [0, 1]
        # option for right segment
        op4 = copy.deepcopy(opvert)
        op4.start_point = [0, 1]
        op4.end_point = [0, 0]
        # building boundary grid
        bgrid = hmscript.BuildBoundaryGrid([op1, op2, op3, op4])

    """
    inp = []
    if not isinstance(opts, list):
        opts2 = [opts]
    else:
        opts2 = opts
    for op in opts2:
        d = {}
        d['source'] = op.contour_id
        d['partition'] = op.partition
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
        else:
            raise ValueError("Invalid bnd_stepping: " + str(op.bnd_stepping))
        d['mesh_cont_step'] = op.bnd_step
        d['algo_acute'] = op.range_angles[0]
        d['algo_right'] = op.range_angles[1]
        d['algo_straight'] = op.range_angles[2]
        d['algo_reentr'] = op.range_angles[3]
        if (op.start_point is not None and op.end_point is not None):
            d['start'] = Point2(*op.start_point)
            d['end'] = Point2(*op.end_point)
        d['force_conf'] = op.force_conformal
        inp.append(d)

    c = com.gridcom.BuildBoundaryGrid({"opt": inp})
    try:
        flow.exec_command(c)
        return c._get_added_names()[0][0]
    except Exception:
        raise ExecError('build_boundary_grid')
