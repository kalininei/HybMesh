"2d grids operations and modifications"
from hybmeshpack import com
from hybmeshpack.hmscript import flow, hmscriptfun
import o2info
from datachecks import (icheck, List, UListOr1, Bool, UList, Point2D,
                        Grid2D, OneOf, Float, Tuple, ACont2D, CompoundList)


@hmscriptfun
def exclude_contours(grid, conts, what):
    """Builds a grid by excluding contour area from existing grid.

    :param grid: source grid identifier

    :param conts: contour or list of contour/grid identifiers for exclusion.

    :param str what: ``"inner"``/``"outer"``.
       Defines what part of ``conts`` domain should be excluded.

    :returns: new grid identifier

    See details in :ref:`substract-area`.

    .. note::

       All contours from **cont** list are excluded consecutively.
       If you want to exclude multiply connected domain area you should first
       assemble multiply connected domain from the list of singly
       connected ones using :func:`unite_contours` procedure.

    Example:
       .. literalinclude:: ../../testing/py/fromdoc/ex_exclude.py
           :start-after: START OF EXAMPLE
           :end-before: END OF EXAMPLE

    """
    icheck(0, Grid2D())
    icheck(1, UListOr1(ACont2D()))
    icheck(2, OneOf("inner", "outer"))

    if not isinstance(conts, list):
        conts = [conts]
    if what == "inner":
        exa = True
    elif what == "outer":
        exa = False
    c = com.gridcom.ExcludeContours({"grid_name": grid,
                                     "cont_names": conts,
                                     "is_inner": exa})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def unite_grids(base_grid, over_grids, empty_holes=False, fix_bnd=False,
                zero_angle_approx=0, buffer_fill='3'):
    """Makes grids superposition.

    :param base_grid: basic grid identifier

    :param list-of-tuples over_grids: sequence of grids for superposition as
      ``[(grid_id, buffer), () ...]`` where
      ``grid_id`` is an superposed grid identifier,
      ``buffer`` - size of the buffer for current union

    :param bool empty_holes: keep all empty zones
      (in case of multiple connectivity)
      of imposed grids in the resulting grid.

    :param bool fix_bnd: whether to fix all boundary nodes

    :param positive-degree zero_angle_approx:
      defines deviation from the straight angle which is considered
      insignificant. Grid boundary vertices which provide insignificant
      contour turns could be moved in order to obtain better result.
      Makes sense only if ``fix_bnd = False``.

    :param str buffer_fill: type of grid in a buffer.

      * ``"3"`` - triangle grid
      * ``"4"`` - mostly quadrangle grid

    :return: identifier of the newly created grid

    See detailed options description in :ref:`gridimp`.

    Example:
       .. literalinclude:: ../../testing/py/fromdoc/ex_unite_grids.py
           :start-after: START OF EXAMPLE
           :end-before: END OF EXAMPLE

    """
    icheck(0, Grid2D())
    icheck(1, UList(Tuple(Grid2D(), Float())))
    icheck(2, Bool())
    icheck(3, Bool())
    icheck(4, Float(within=[-1., 180., '[]']))
    icheck(5, OneOf('3', '4'))

    args = {"base": base_grid, "empty_holes": empty_holes,
            "angle0": zero_angle_approx,
            "fix_bnd": fix_bnd, "plus": [],
            "filler": buffer_fill}
    for ig in over_grids:
        args["plus"].append({"name": ig[0], "buf": ig[1], "den": 7})
    c = com.gridcom.UniteGrids(args)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def unite_grids1(base_grid, secondary_grid, buffer_size,
                 empty_holes=False, fix_bnd=False,
                 zero_angle_approx=0, buffer_fill='3'):
    """Makes superposition of two grids.

    :param base_grid: basic grid identifier

    :param secondary_grid: superposed grid identifier

    :param float buffer_size: size of the buffer zone

    :param bool empty_holes: keep all empty zones
      (in case of multiple connectivity)
      of imposed grids in the resulting grid.

    :param bool fix_bnd: whether to fix all boundary nodes

    :param positive-degree zero_angle_approx:
      defines deviation from the straight angle which is considered
      insignificant. Grid boundary vertices which provide insignificant
      contour turns could be moved in order to obtain better result.
      Makes sense only if ``fix_bnd = False``.

    :param str buffer_fill: type of grid in a buffer.

      * ``"3"`` - triangle grid
      * ``"4"`` - mostly quadrangle grid

    :return: identifier of the newly created grid

    Same as :func:`unite_grids` without chain superposition support.
    """
    icheck(0, Grid2D())
    icheck(1, Grid2D())
    icheck(2, Float())
    icheck(3, Bool())
    icheck(4, Bool())
    icheck(5, Float(within=[-1., 180., '[]']))
    icheck(6, OneOf('3', '4'))

    return unite_grids(base_grid, [(secondary_grid, buffer_size)],
                       empty_holes, fix_bnd, zero_angle_approx,
                       buffer_fill)


@hmscriptfun
def remove_cells(grid, cont, what):
    """ Removes grid cells by a bounding contour.

    :param grid: basic grid identifier

    :param cont: bounding area, grid or contour identifier

    :param str what: remove cells which are

      * ``"fully_inside"``,
      * ``"fully_outside"``,
      * ``"partly_inside"``,
      * ``"partly_outside"``,
      * ``"cross"``, -- crossed by **cont** edges,
      * ``"no_cross"``, -- not crossed by **cont** edges.

    :return: identifier of the newly created grid

    See details in :ref:`remove-cells`.

    """
    icheck(0, Grid2D())
    icheck(1, ACont2D())
    icheck(2, OneOf('fully_inside', 'fully_outside',
        'partly_inside', 'partly_outside', 'cross', 'no_cross'))

    args = {"grid_name": grid, "cont_name": cont, "what": what}
    c = com.gridcom.RemoveCells(args)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def inscribe_grid(grid, cont, where, buffer_size,
        keep_cont=True, zero_angle=0, buffer_fill='3'):
    """ Inscribes grid into an area.

    :param grid: input grid identifier

    :param cont: bounding area, grid or contour identifier

    :param str where: inscribe inside/outside **cont** area.

      * ``"inside"``,
      * ``"outside"``.

    :param float buffer_size: offset distance from **cont**
       area which defines a buffer zone.

    :param bool keep_cont: whether use **cont** segmentation as is or
      repart it.

    :param degree zero_angle:
      defines deviation from the straight angle which is considered
      insignificant while 1D resegmentation. -1 makes all existed vertices
      significent; 180 -- allows complete resegmentation.

    :param str buffer_fill: one of

      * ``"3"`` - triangulate buffer,
      * ``"4"`` - build recombined grid,
      * ``"no"`` - do not fill the buffer with unstructured grid.

    :return: identifier of the newly created grid


    See details in :ref:`inscribe-grid`.

    """
    icheck(0, Grid2D())
    icheck(1, ACont2D())
    icheck(2, OneOf('inside', 'outside'))
    icheck(3, Float(grthan=0.0))
    icheck(4, Bool())
    icheck(5, Float(within=[-1.0, 180., '[]']))
    icheck(6, OneOf('3', '4', 'no'))

    args = {"grid_name": grid, "cont_name": cont, "where": where,
        "buffer": buffer_size, "keep_cont": keep_cont,
        "zero_angle": zero_angle, "buffer_fill": buffer_fill}
    c = com.gridcom.InscribeGrid(args)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def insert_grid_constraints(grid, conts=[], sites=[], buffer_size=1.,
        keep_cont=True, zero_angle=0, buffer_fill='3'):
    """ Inserts a line or site constraint into a grid.

    :param grid: input grid identifier

    :param conts: contour (or grid) identifier or list of identifiers 
       for line constraints

    :param sites: list of site constraints given as
       ``[size0, [x0, y0], size1, [x1, y1], ...]`` where *size0,1,...* is the
       desired edge size in a given site. Set it to -1 or ``'auto'``
       for automatic size detection.

    :param float buffer_size: offset distance from **conts** and **sites**
       which defines a buffer zone.

    :param keep_cont: whether line constraints should be used as they are or
       remeshed with a proper step.

    :param degree zero_angle:
      defines deviation from the straight angle which is considered
      insignificant while 1D resegmentation. -1 makes all existed vertices
      significant.

    :param str buffer_fill: one of

      * ``"3"`` - triangulate buffer,
      * ``"4"`` - build recombined grid,
      * ``"no"`` - do not fill the buffer with unstructured grid.

    :return: identifier of the newly created grid

    See details in :ref:`insert-constraints`.

    """
    i = 0;
    while i<len(sites):
        if sites[i] == 'auto':
            sites[i] = -1
        i += 2
    if not isinstance(conts, list):
        conts = [conts]

    icheck(0, Grid2D())
    icheck(1, UList(ACont2D()))
    icheck(2, CompoundList(Float(), Point2D()))
    icheck(3, Float(grthan=0.0))
    icheck(4, Bool())
    icheck(5, Float(within=[-1.0, 180., '[]']))
    icheck(6, OneOf('3', '4', 'no'))

    args = {"grid_name": grid, "conts": conts, "sites": [], 
        "buffer": buffer_size, "keep_cont": keep_cont,
        "zero_angle": zero_angle, "buffer_fill": buffer_fill}
    it = iter(sites)
    for dist, pnt in zip(it, it):
        args['sites'].append(dist)
        args['sites'].append(pnt[0])
        args['sites'].append(pnt[1])
    c = com.gridcom.InsertConstraints(args)
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def map_grid(base_grid, target_contour, base_points, target_points,
             snap="no", project_to="line", btypes="from_grid",
             algo="inverse_laplace",
             is_reversed=False, return_invalid=False):
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

    :param bool is_reversed: shows whether target contour
       should be treated in reversed order while building boundary mapping.

    :param bool return_invalid: if this flag is on
       then the procedure will return a grid even if it is not valid
       (has self-intersections). Such grids could be exported to
       simple formats (like vtk or tecplot) in order to detect
       bad regions and give user a hint of how to adopt
       input data to gain an acceptable result.

       .. warning:: Never use invalid grids for further operations.

    :returns: identifier of newly created grid

    """
    icheck(0, Grid2D())
    icheck(1, ACont2D())
    icheck(2, List(Point2D(), minlen=1))
    icheck(3, List(Point2D(), llen=len(base_points)))
    icheck(4, OneOf('no', 'add_vertices', 'shift_vertices'))
    icheck(5, OneOf('line', 'vertex', 'corner'))
    icheck(6, OneOf('from_grid', 'from_contour'))
    icheck(7, OneOf('inverse_laplace', 'direct_laplace'))
    icheck(8, Bool())
    icheck(9, Bool())

    # project_to option treatment
    if project_to == "line":
        pass
    elif project_to == "vertex":
        base_points = [o2info.get_point(base_grid, vclosest=p,
                                        only_contour=True)
                       for p in base_points]
        target_points = [o2info.get_point(target_contour, vclosest=p,
                                          only_contour=True)
                         for p in target_points]
    elif project_to == "corner":
        base_points = [o2info.get_point(base_grid, cclosest=p,
                                        only_contour=True)
                       for p in base_points]
        target_points = [o2info.get_point(target_contour, cclosest=p,
                                          only_contour=True)
                         for p in target_points]
    c = com.gridcom.MapGrid({"base": base_grid,
                             "target": target_contour,
                             "base_points": base_points,
                             "target_points": target_points,
                             "snap": snap,
                             "algo": algo,
                             "btypes": btypes,
                             "is_reversed": is_reversed,
                             "return_invalid": return_invalid})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def snap_grid_to_contour(gid, cid, gstart, gend, cstart, cend, algo="add"):
    """
    Snaps grid boundary subsection to contour subsection.
    See details in :ref:`snapgrid`.

    :param gid: source grid identifier.

    :param cid: target contour identifier.

    :param gstart:

    :param gend: grid boundary start/end points as ``[x, y]``.

    :param cstart:

    :param cend: contour start/end points as ``[x, y]``.

    :param str algo: Intermediate grid vertices will be:

        * ``'add'`` - projected to contour;
        * ``'shift'`` - shifted to closest contour vertices.

    :returns: new grid identifier.

    Given end points will be projected to closest objects vertices.
    To define the whole closed contour let start point
    be equal to respective end point.

    """
    icheck(0, Grid2D())
    icheck(1, ACont2D())
    icheck(2, Point2D())
    icheck(3, Point2D())
    icheck(4, Point2D())
    icheck(5, Point2D())
    icheck(6, OneOf("add", "shift"))
    c = com.gridcom.SnapToContour({"grid": gid, "cont": cid,
                                   "gp1": gstart, "gp2": gend,
                                   "cp1": cstart, "cp2": cend,
                                   "algo": algo})
    flow.exec_command(c)
    return c.added_grids2()[0]


@hmscriptfun
def heal_grid(gid, simplify_boundary=30, convex_cells=-1):
    """ Set of procedures for simplification of grid geometry

    :param gid: identifier or (list of identifiers) of the grid

    :param float simplify_boundary: angle (deg) in [0, 180].

    :param float convex_cells: angle (deg) in [0, 360]. Maximum
      allowed angle in a non-triangle cell.

    :return: None

    If **simplify_boundary** parameter is non-negative then edges which

    * are boundary,
    * belong to the same grid cell,
    * form an angle no more then ``simplify_boundary`` degree,

    will be merged if possible. If ``simplify_boundary=0`` then only edges
    lying on the same line are considered, ``simplify_boundary=180``
    leads to merging of all doubled boundary edges,
    ``simplify_boundary=-1`` ignores this simplification option.

    if **convex_cells** is non negative than all cells which provide
    angles greater or equal to this value will be splitted.
    If ``convex_cells=180`` then all concave cells includings those
    whith hanging nodes will be processed.
    If you set ``convex_cells=0`` then all cells will be splitted
    up to triangles.

    Subprocedures order is:

    * boundary simplification
    * remove concave cells

    """
    icheck(0, UListOr1(Grid2D()))
    icheck(1, Float(within=[-1., 180., '[]']))
    icheck(2, Float(within=[-1., 360., '[]']))

    if isinstance(gid, list):
        for g in gid:
            heal_grid(g, simplify_boundary, convex_cells)
        return
    c = com.gridcom.HealGrid({"name": gid,
                              "simp_bnd": simplify_boundary,
                              "convex": convex_cells})
    flow.exec_command(c)
