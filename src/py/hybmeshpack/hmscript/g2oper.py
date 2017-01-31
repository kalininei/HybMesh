"2d grids operations and modifications"
from hybmeshpack import com
from hybmeshpack.hmscript import flow, ExecError
import o2info


def exclude_contours(grid, conts, what):
    """Builds a grid by excluding contour area from existing grid.

    :param grid: source grid identifier

    :param conts: contour or list of contours/grids identifiers for exclusion.

    :param str what: ``"inner"``/``"outer"``.
       Describes what part of ``conts`` domain exclude

    :returns: new grid identifier

    :raises: ValueError, ExecError

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
    if not isinstance(conts, list):
        conts = [conts]
    if what == "inner":
        exa = True
    elif what == "outer":
        exa = False
    else:
        raise ValueError("Invalid 'what' parameter: %s" % str(what))
    c = com.gridcom.ExcludeContours({"grid_name": grid,
                                     "cont_names": conts,
                                     "is_inner": exa})
    try:
        flow.exec_command(c)
        return c.added_grids2[0]
    except:
        raise ExecError('exclude_contours')


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

    :raises: hmscript.ExecError

    See detailed options description in :ref:`gridimp`.

    Example:
       .. literalinclude:: ../../testing/py/fromdoc/ex_unite_grids.py
           :start-after: START OF EXAMPLE
           :end-before: END OF EXAMPLE

    """
    args = {"base": base_grid, "empty_holes": empty_holes,
            "angle0": zero_angle_approx,
            "fix_bnd": fix_bnd, "plus": [],
            "filler": buffer_fill}
    for ig in over_grids:
        args["plus"].append({"name": ig[0], "buf": ig[1], "den": 7})
    if buffer_fill not in ['3', '4']:
        raise ValueError("Invalid buffer_fill option")
    c = com.gridcom.UniteGrids(args)
    try:
        flow.exec_command(c)
        return c.added_grids2()[0]
    except Exception:
        raise ExecError('unite_grids')


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

    :raises: ValueError, hmscript.ExecError

    :returns: identifier of newly created grid

    """
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
    try:
        flow.exec_command(c)
        return c.added_grids2[0]
    except:
        raise ExecError("map_grid")


def heal_grid(gid, simplify_boundary=30, convex_cells=-1):
    """ Set of procedures for simplification of grid geometry

    :param gid: identifier or (list of identifiers) of the grid

    :param float simplify_boundary: angle (deg) in [0, 180].

    :param float convex_cells: angle (deg) in [0, 180]

    :return: None

    :raises: ValueError, ExecError

    If **simplify_boundary** parameter is non-negative then edges which

    * are boundary,
    * belong to the same grid cell,
    * form an angle no more then ``simplify_boundary`` degree,

    will be merged if possible. If ``simplify_boundary=0`` then only edges
    lying on the same line are considered, ``simplify_boundary=180``
    leads to merging of all doubled boundary edges,
    ``simplify_boundary=-1`` ignores this simplification option.

    if ``convex_cells`` is non negative than all concave cells
    will be turned into convex ones by their division. Paramter represents
    concave angle at which procedure will be executed. ``0`` provides
    elimination of all concave segments. ``180`` - only degenerate ones.

    Subprocedures order is:

    * boundary simplification
    * remove concave cells

    """
    if isinstance(gid, list):
        for g in gid:
            heal_grid(g, simplify_boundary, convex_cells)
        return
    c = com.gridcom.HealGrid({"name": gid,
                              "simp_bnd": simplify_boundary,
                              "convex": convex_cells})
    try:
        flow.exec_command(c)
    except:
        raise ExecError('heal_grid')