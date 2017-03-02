"2D grids construction procedures"
import copy
from hybmeshpack import com
from hybmeshpack.hmscript import flow, hmscriptfun
import o2info
from datachecks import (icheck, List, UListOr1, Bool, Point2D, ListOr1, ZType,
                        UInt, ACont2D, OneOf, Float, NoneOr, InvalidArgument,
                        Cont2D, Or, IncList, CompoundList, Any)


# Prototype grids
@hmscriptfun
def add_unf_rect_grid(p0=[0, 0], p1=[1, 1], nx=3, ny=3,
                      custom_x=[], custom_y=[], bnd=0):
    """Builds rectangular grid.

    :param list-of-floats p0:

    :param list-of-floats p1: bottom left, top right points in [x, y] format.

    :param int nx:

    :param int ny: partition in x and y directions.

    :param float-or-list-of-floats custom_x:

    :param float-or-list-of-floats custom_y: custom x and y coordinates

    :param int-or-list-of-int: boundary types for bottom, right, top, left
       rectangle sides

    :returns: created grid identifier

    Builds a grid in a rectangular area formed by points **p0** and **p1**.
    **nx** and **ny** provide grid partition in x and y direction.

    If **custom_x**/**custom_y** is given by a single float value
    than it shows a step size in respective direction,
    hence values given by **nx**/**ny** parameters will be omitted.

    If **custom_x**/**custom_y** is given by a list of increasing floats
    it explicitly shows the partition in respective direction.
    In the latter case the respective **p0**, **p1** coordinates
    will also be ignored.

    Use :func:`partition_segment` to conveniently define **custom\_** fields
    if needed.
    """
    icheck(0, Point2D())
    icheck(1, Point2D(grthan=p0))
    icheck(2, UInt(minv=1))
    icheck(3, UInt(minv=1))
    icheck(4, Or(Float(grthan=0.), IncList(Float())))
    icheck(5, Or(Float(grthan=0.), IncList(Float())))
    icheck(6, ListOr1(ZType(), llen=4))

    bnd = bnd[:4] if isinstance(bnd, list) else [bnd, bnd, bnd, bnd]
    custom_x = custom_x if isinstance(custom_x, list) else [custom_x]
    custom_y = custom_y if isinstance(custom_y, list) else [custom_y]
    c = com.gridcom.AddUnfRectGrid({"p0": p0, "p1": p1,
                                    "nx": nx, "ny": ny,
                                    "custom_x": custom_x,
                                    "custom_y": custom_y,
                                    "bnds": bnd})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def add_unf_circ_grid(p0, rad=1.0, na=8, nr=4, coef=1.0, is_trian=True,
                      custom_rads=[], custom_archs=[], bnd=0):
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

    :param float-or-list-of-floats custom_rads:
        user defined radious partition

    :param float-or-list-of-floats custom_archs:
        user defined arch partition

    :param int bnd: boundary type for outer contour

    :returns: created grid identifier

    Creates a radial grid with the center in **p0**.

    If **custom_rads** is given as a single value it will be used
    as a constant step along radial axis hence **nr** and **coef** arguments
    will be ignored. If it is given as a list of increasing values
    starting from zero
    then it is parsed as explicit radius partition. Hence the
    last entry of this list will be the radius of the resulting grid
    and **rad** parameter will also be ignored.

    If **custom_archs** is given as a single value it shows the
    constant step along outer arch and **na** will be ignored.
    If it is an increasing list of floats, it shows partition of
    outer arch. It can be given in degrees or radians or any other
    relative units. Program treats **custom_archs[-1]**-**custom_archs[0]**
    difference as a full circle length and normalizes all other entries
    to this length. First and last entries of this array provides
    the same arch segment (last = first + 2*pi) hence to
    get partition of n segments you should define n+1 entries.

    Use :func:`partition_segment` to conveniently define **custom\_** fields
    if needed.
    """
    icheck(0, Point2D())
    icheck(1, Float(grthan=0.0))
    icheck(2, UInt(minv=3))
    icheck(3, UInt(minv=1))
    icheck(4, Float(grthan=0.0))
    icheck(5, Bool())
    icheck(6, Or(Float(grthan=0.0), IncList(Float())))
    icheck(7, Or(Float(grthan=0.0), IncList(Float())))
    icheck(8, ZType())

    custom_rads = custom_rads if isinstance(custom_rads, list)\
        else [custom_rads]
    custom_archs = custom_archs if isinstance(custom_archs, list)\
        else [custom_archs]
    c = com.gridcom.AddUnfCircGrid({
        "p0": p0, "rad": rad,
        "na": na, "nr": nr,
        "coef": coef,
        "is_trian": is_trian,
        "custom_r": custom_rads,
        "custom_a": custom_archs,
        "bnd": bnd})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def add_unf_ring_grid(p0, radinner, radouter,
                      na, nr, coef=1.0, bnd=0):
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

    :param int-or-list-of-int bnd: boundary types for inner and outer
       ring boundaries

    :return: created grid identifier

    """
    if radinner > radouter:
        radinner, radouter = radouter, radinner
    icheck(0, Point2D())
    icheck(1, Float(grthan=0.0))
    icheck(2, Float(grthan=0.0))
    icheck(3, UInt(minv=3))
    icheck(4, UInt(minv=1))
    icheck(5, Float(grthan=0.0))
    icheck(6, ListOr1(ZType(), llen=2))

    bnd = bnd[:2] if isinstance(bnd, list) else [bnd, bnd]
    c = com.gridcom.AddUnfRingGrid({
        "p0": p0,
        "radinner": radinner, "radouter": radouter,
        "na": na, "nr": nr,
        "coef": coef,
        "bnd": bnd})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def add_unf_hex_grid(area, cell_radius, strict=False):
    """ Builds grid with regular hexagonal cells

    :param area: defines meshing area. If given as ``[[x, y], radius]``
       then represents hexagonal area; if ``[[x0, y0], [x1, y1]]`` then
       it is a rectangle defined by bottom-left and top-right points.

    :param float cell_radius: radius of hexagonal cell.

    :param bool strict: forces grid stretch
      to guarantee that all outer rectangle corners lie in the
      centers of cells.

    See details in :ref:`hexgrid`
    """
    icheck(0, List(Any()))
    icheck(0, CompoundList(Point2D(),
                           Or(Float(grthan=0.0),
                              Point2D(grthan=area[0])),
                           llen=2))
    icheck(1, Float(grthan=0.))
    icheck(2, Bool())

    simpar = [area[0][0], area[0][1]]
    if isinstance(area[1], list):
        simpar.append(area[1][0])
        simpar.append(area[1][1])
    else:
        simpar.append(area[1])
    args = {"area": simpar, "crad": cell_radius, "strict": strict}
    c = com.gridcom.AddUnfHexGrid(args)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def add_triangle_grid(p0, p1, p2, nedge, bnd=0):
    """Creates structured grid in triangle area

    :param list-of-floats p0:

    :param list-of-floats p1:

    :param list-of-floats p2: triangle vertices in [x, y] format

    :param int nedge: partition of triangle edges

    :param int-or-list-of-int bnd: boundary types for outer contour

    :return: identifier of newly created grid

    Resulting grid will contain quadrangle cells everywhere except
    area near ``p0``-``p2`` edge where triangle cells will be built.
    """
    icheck(0, Point2D())
    icheck(1, Point2D(noteq=[p0]))
    icheck(2, Point2D(noteq=[p0, p1]))
    icheck(3, UInt(minv=1))
    icheck(4, ListOr1(ZType(), llen=3))

    bnd = bnd[:3] if isinstance(bnd, list) else [bnd, bnd, bnd]
    c = com.gridcom.AddTriGrid({
        "vertices": [p0, p1, p2], "nedge": nedge, "bnd": bnd})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def add_custom_rect_grid(algo, left, bottom, right=None, top=None,
                         hermite_tfi_w=[1.0, 1.0, 1.0, 1.0],
                         return_invalid=False):
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

    :param bool return_invalid: if this flag is on
       then the procedure will return a grid even if it is not valid
       (has self-intersections). Such grids could be exported to
       simple formats (like vtk or tecplot) in order to detect
       bad regions and give user a hint of how to adopt
       input data to gain an acceptable result.

       .. warning:: Never use invalid grids for further operations.

    :return: new grid identifier

    """
    icheck(0, OneOf('linear', 'linear_tfi', 'hermite_tfi', 'inverse_laplace',
                    'direct_laplace', 'orthogonal'))
    icheck(1, Cont2D())
    icheck(2, Cont2D())
    icheck(3, NoneOr(Cont2D()))
    icheck(4, NoneOr(Cont2D()))
    icheck(5, List(Float(), llen=4))
    icheck(6, Bool())

    # call
    args = {'algo': algo,
            'left': left,
            'right': right,
            'bot': bottom,
            'top': top,
            'her_w': hermite_tfi_w,
            'return_invalid': return_invalid}
    c = com.gridcom.AddCustomRectGrid(args)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
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

    """
    icheck(0, Point2D())
    icheck(1, Float(grthan=0.0))
    icheck(2, Float(grthan=0.0))
    icheck(3, Float(within=[0., 1.4, '[)']))
    icheck(4, Float(grthan=0.0))
    icheck(5, OneOf('linear', 'laplace', 'orthogonal_circ', 'orthogonal_rect'))

    # call
    args = {'algo': algo,
            'p0': p0,
            'rad': rad,
            'step': step,
            'sqrside': sqrside,
            'rcoef': rcoef}
    c = com.gridcom.AddCirc4Grid(args)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def stripe(cont, partition, tip='no', bnd=0):
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

    Horizontal partition is taken from contour partition.
    Vertical partition is given by user with ``partition`` list parameter.
    If it starts with non zero value then grid will not contain
    contour nodes as its vertices.

    Use :func:`partition_segment` to define non-equidistant
    **partition** with any desired refinement if needed.
    """
    icheck(0, ACont2D())
    icheck(1, IncList(Float(grthan=0.0)))
    icheck(2, OneOf('no', 'radial'))
    icheck(3, ListOr1(ZType(), llen=4))

    bnd = bnd[:4] if isinstance(bnd, list) else [bnd, bnd, bnd, bnd]
    arg = {"source": cont, "partition": partition, "tip": tip, "bnd": bnd}
    c = com.gridcom.StripeGrid(arg)
    flow.exec_command(c)
    return c.added_grids2()[0]


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
    if (tp == '4'):
        c = com.gridcom.QuadrangulateArea(arg)
    elif (tp == 'pebi'):
        c = com.gridcom.PebiFill(arg)
    else:
        c = com.gridcom.TriangulateArea(arg)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def triangulate_domain(domain, constr=[], pts=[], fill='3'):
    """Builds constrained triangulation within given domain

    :param domain: single or list of closed contours
        representing bounding domain

    :param constr: single or list of contours representing
        triangulation constraints

    :param pts: set of points in ``[len0, [x0, y0], ...]``
        format where ``x, y`` are coordinates of internal vertices
        which should be embedded into the resulting grid,
        ``len`` - size of adjacent cells
    :param str fill: if '3' then triangulates area; '4' runs
        recombination algorithm to make mostly quadrangular mesh

    :return: grid identifier

    A contour tree will be built using all closed contours
    passed as **domain** parameter. Only the interior parts
    of this tree will be meshed. Contours passed by **domain**
    should not intersect each other, but could intersect **constr**
    contours.
    **constr** could contain any set of closed and open contours.

    See details in :ref:`unstructured-meshing`.
    """
    icheck(0, UListOr1(ACont2D()))
    icheck(1, UListOr1(ACont2D()))
    icheck(2, CompoundList(Float(grthan=0.0), Point2D()))
    icheck(3, OneOf('3', '4'))

    if fill == '3':
        return _triquad(domain, constr, pts, '3')
    elif fill == '4':
        return _triquad(domain, constr, pts, '4')


@hmscriptfun
def pebi_fill(domain, constr=[], pts=[]):
    """Builds perpendicular bisector cells in given domain.

    :param domain: single or list of closed contours
        representing bounding domain

    :param constr: single or list of contours representing
        meshing constraints

    :param pts: set of points in ``[len0, [x0, y0], ...]``
        format where ``x, y`` are coordinates of internal vertices
        which should be embedded into the resulting grid,
        ``len`` - size of adjacent cells

    :return: grid identifier

    A contour tree will be built using all closed contours
    passed as **domain** parameter. Only the interior parts
    of this tree will be meshed. Contours passed by **domain**
    should not intersect each other, but could intersect **constr**
    contours.

    Routine can produce concave cells (as a result of bad size
    control or near the concave domain boundary vertices).
    Use :func:`heal_grid` routine with ``convex_cells`` option to fix this.

    See details in :ref:`unstructured-meshing`.
    """
    icheck(0, UListOr1(ACont2D()))
    icheck(1, UListOr1(ACont2D()))
    icheck(2, CompoundList(Float(grthan=0.0), Point2D()))
    return _triquad(domain, constr, pts, 'pebi')


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

    Class contains functions to define uniform and increasing
    vertical partitions. However user could give his own custom partition by
    explicit assigning of **partition** attribute.
    Use :func:`partition_segment` for advanced control of vertical
    partition.
    """

    @hmscriptfun
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

    @hmscriptfun
    def uniform_partition(self, fullh, n):
        """Sets uniform boundary grid vertical `partition` with
        the full depth of ``fullh`` and ``n`` total layers
        """
        h = fullh / n
        self.partition = [i * h for i in range(n + 1)]

    @hmscriptfun
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


@hmscriptfun
def build_boundary_grid(opts):
    """Builds a boundary grid near contour

    :params opts: single or list of :class:`BoundaryGridOptions` objects

    :returns: identifier of the newly created grid.

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
        if op.direction == "left":
            d['direction'] = 1
        elif op.direction == "right":
            d['direction'] = -1
        else:
            raise InvalidArgument("Invalid direction: " + str(op.direction))
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
            raise InvalidArgument(
                "Invalid bnd_stepping: " + str(op.bnd_stepping))

        if op.bnd_stepping != "incremental":
            if isinstance(op.bnd_step, list):
                raise InvalidArgument("list values for bnd_step are available "
                                      "only for incremental stepping")
            else:
                d['step_start'], d['step_end'] = 1., 1.
                d['mesh_cont_step'] = op.bnd_step
        else:
            d['mesh_cont_step'] = 0.0
            if not isinstance(op.bnd_step, list) or len(op.bnd_step) < 2:
                raise InvalidArgument(
                    "Incremental stepping requires list[2] as bnd_step")
            d['step_start'] = op.bnd_step[0]
            d['step_end'] = op.bnd_step[1]

        d['algo_acute'] = op.range_angles[0]
        d['algo_right'] = op.range_angles[1]
        d['algo_straight'] = op.range_angles[2]
        d['algo_reentr'] = op.range_angles[3]
        d['start'], d['end'] = None, None
        if (op.start_point is not None and op.end_point is not None):
            d['start'], d['end'] = op.start_point, op.end_point
            if op.project_to == "line":
                pass
            elif op.project_to == "vertex":
                d['start'] = o2info.get_point(
                    op.contour_id, vclosest=d['start'], only_contour=True)
                d['end'] = o2info.get_point(
                    op.contour_id, vclosest=d['end'], only_contour=True)
            elif op.project_to == "corner":
                d['start'] = o2info.get_point(
                    op.contour_id, cclosest=d['start'], only_contour=True)
                d['end'] = o2info.get_point(
                    op.contour_id, cclosest=d['end'], only_contour=True)
            else:
                raise InvalidArgument(
                    "Unknown `project_to` = %s" % op.project_to)

        d['force_conf'] = op.force_conformal
        # check partition
        for i in range(len(op.partition) - 1):
            if (op.partition[i] >= op.partition[i + 1]):
                op.partition = [0]
                break
        if len(op.partition) < 2:
            raise InvalidArgument("Invalid partition")
        d['partition'] = op.partition
        inp.append(d)

    c = com.gridcom.BuildBoundaryGrid({"opt": inp})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def simple_boundary_grid(cont, partition, direction, pstart=None, pend=None,
                         range_angles=[40, 125, 235, 275]):
    """Builds a singly-connected boundary grid near contour

    :params cont: source contour (or grid) identifier

    :params list-of-float partition: partition in perpendicular direction.

    :params str direction: 'left'/'right'

    :params pstart:

    :params pend: points in [x, y] format which define
      the exact segment of the contour for building grid.
      If both are None hence whole contour (or all subcontours) will be used.

    :params range_angles: list of 4 angle values (deg) which define algorithms
      for contour bends treatment.

    :returns: identifier of the newly created grid.

    This is a wrapper for a :func:build_boundary_grid with simplified
    interface. It allows to build a boundary grid with constant partition
    options using existing contour segmentation for horizontal stepping.
    """
    icheck(0, ACont2D())
    icheck(1, IncList(Float(), startfrom=0.0))
    icheck(2, OneOf('left', 'right'))
    icheck(3, NoneOr(Point2D()))
    icheck(4, NoneOr(Point2D(noteq=pstart)))
    icheck(5, IncList(Float(within=[0, 360, '[]'])))

    bo = BoundaryGridOptions(
        cont, partition, direction, bnd_stepping='no',
        range_angles=range_angles, start_point=pstart,
        end_point=pend)
    return build_boundary_grid([bo])
