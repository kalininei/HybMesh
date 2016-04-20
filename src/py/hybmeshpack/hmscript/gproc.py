import copy
from hybmeshpack import com
from hybmeshpack.hmscript import flow, data
from hybmeshpack.basic.geom import Point2
from . import ExecError


def remove_geom(objs):
    """ Completely removes object or list of objects

    Args:
       objs: identifier or list of identifiers of removing objects

    """
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.RemoveGeom({"names": ob})
    flow.exec_command(c)


def remove_all():
    """ Completely removes all grids, contours and boundary types
    """
    flow.to_zero_state()


def remove_all_but(objs):
    """ Removes all geometry objects except for listed ones

    Args:
       objs: identifier or list of identifiers of objects
       which should not be removed

    """
    if not isinstance(objs, list):
        objs = [objs]

    all_obj = data.get_grid_names()
    all_obj.extend(data.get_ucontour_names())
    all_obj.extend(data.get_grid3_names())
    all_obj = [x for x in all_obj if x not in objs]
    remove_geom(all_obj)


def move_geom(objs, dx, dy):
    """ Moves a list of objects

    Args:
       objs: identifier or list of identifiers of moving 2d objects

       dx, dy (float): shifts in x and y direction

    """
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.MoveGeom({"names": ob, "dy": dy, "dx": dx})
    flow.exec_command(c)


def rotate_geom(objs, angle, pc=[0.0, 0.0]):
    """ Rotates group of objects

    Args:
       objs: identifier or list of identifiers of rotating 2d objects

       angle (float): degree of rotation. Positive angle corresponds to
       counterclockwise rotation

    Kwargs
       pc (list-of-float): center of rotation
    """
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.RotateGeom({"names": ob, "angle": angle,
                               "p0": Point2(*pc)})
    flow.exec_command(c)


def scale_geom(objs, xpc, ypc, refp=[0.0, 0.0]):
    """ Scales objects

    Args:
       objs: identifier or list of identifiers of scaling 2d objects

       xpc, ypc (float): percentages of scaling in x and y directions

       refp (list-of-float): reference point which stays
       fixed after transformation
    """
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.ScaleGeom({"names": ob, "xpc": xpc, "ypc": ypc,
                              "p0": Point2(*refp)})
    flow.exec_command(c)


def copy_geom(objs):
    """ Creates deep copies of objects

    Args:
       objs: identifier or list of identifiers of objects to copy

    Returns:
       list of identifiers of copied objects
    """
    if isinstance(objs, list):
        ret = []
        for s in objs:
            ret.append(copy_geom(s)[0])
        return ret

    c = com.objcom.CopyGeom({"names": [objs], "newnames": [objs]})
    flow.exec_command(c)
    if len(c._get_added_names()[0]) > 0:
        return c._get_added_names()[0]
    else:
        return c._get_added_names()[1]


def reflect_geom(objs, pnt1, pnt2):
    """ Makes a reflection of geometry objects
    over a  given line.

    Args:
       objs: identifier or list of identifiers of 2d objects to reflect

       pnt1, pnt2: points in [x, y] format which define a line to reflect over

    """
    if not isinstance(objs, list):
        ob = [objs]
    else:
        ob = objs
    p1, p2 = Point2(*pnt1), Point2(*pnt2)
    c = com.objcom.ReflectGeom({"names": ob, "p1": p1, "p2": p2})
    flow.exec_command(c)


def extrude_grid(obj, zcoords, bottombc=0, topbc=0, sidebc=None):
    """ Creates 3D grid by extrusion of 2D grid along z-axis

    :param obj: 2d grid identifier

    :param list-of-floats zcoords: increasing vector of z values
      which will be used to create 3d points

    :param bottombc:

    :param topbc: values which define boundary features of
      3d grid at ``z=min(zcoords)`` and ``z=max(zcoords)``
      surfaces respectively.
      Could be either a single boundary identifier for a whole
      surface or a function: ``(float x, float y, int cell_index)->bindex``
      which takes central cell point x, y
      coordinates and cell index as arguments and returns boundary type
      (see example below).

    :param sidebc: defines boundary features for side surfaces.

      + If None than boundary types will be taken from corresponding
        edges of 2D grid
      + If single boundary identifier then whole side surface will
        have same boundary type

    :returns: 3D grid identifier

    :raises: ValueError, hmscript.ExecError

    Example:

      .. literalinclude:: ../../testing/py/fromdoc/ex_extrude.py
          :start-after: vvvvvvvvvvvvvvvvvvvvvvvv
          :end-before: ^^^^^^^^^^^^^^^^^^^^^^^^
    """
    # zcoords is strictly increasing vector
    if len(zcoords) < 2:
        raise ValueError("Invalid zcoords vector")
    for i in range(1, len(zcoords)):
            if zcoords[i - 1] >= zcoords[i]:
                raise ValueError("Invalid zcoords vector")
    # calculate boundary types
    if not isinstance(bottombc, int) or not isinstance(topbc, int):
        _, _, g = data.get_grid(name=obj)
        cc_pnt = g.central_points()
        central_coord = [[p.x, p.y, ind] for (ind, p) in enumerate(cc_pnt)]
    if isinstance(bottombc, int):
        bbot = [bottombc]
    elif hasattr(bottombc, '__call__'):
        bbot = [bottombc(x, y, i) for [x, y, i] in central_coord]
    else:
        raise ValueError("Invalid bottombc type")
    if isinstance(topbc, int):
        btop = [topbc]
    elif hasattr(topbc, '__call__'):
        btop = [topbc(x, y, i) for [x, y, i] in central_coord]
    else:
        raise ValueError("Invalid topbc type")
    if sidebc is None:
        bside = None
    elif isinstance(sidebc, int):
        bside = sidebc
    else:
        raise ValueError("Invalid sidebc type")

    c = com.grid3dcom.ExtrudeZ({"base": obj,
                                "zvals": copy.deepcopy(zcoords),
                                "bside": bside,
                                "btop": btop,
                                "bbot": bbot})
    try:
        flow.exec_command(c)
        return c._get_added_names()[2][0]
    except:
        raise ExecError("extrusion")


def revolve_grid(obj, p1, p2, n_phi=None,
                 phi=None, btype1=0, btype2=0, merge_central=False):
    """ Creates 3D grid by revolving of 2D grid around defined vector

    :param obj: 2d grid identifier

    :param p1:

    :param p2: points in [x, y] format which defines vector of rotation

    :param int n_phi: partition along circular coordinate.
       If this parameter is defined then [0, 360] range will be divided
       into equal parts and full revolution solid will be build.

    :param list-of-floats phi: increasing vector defining
       custom partition of angular range.
       This parameter will be processed if **n_phi** is None.
       If the last value of **phi** is not equal to first one
       plus 360 degree than
       incomplete revolution solid will be built.

    :param btype1:

    :param btype2: boundary identifiers for surfaces which will be build
       as a result of incomplete rotation at end values of **phi** vector.

    :param bool merge_central: if rotation vector coincides
       with boundary edges of input grid then this parameter
       defines whether central cells derived from the revolution
       of boundary cells should be merged into one
       complex finite volume (True) or left as they are (False).

    :returns: 3D grid identifier

    :raises: ValueError, hmscript.ExecError

    All points of input grid should lie to the one side of rotation
    vector.

    """
    #check input
    if not isinstance(p1, list) or len(p1) != 2 or\
            not isinstance(p2, list) or len(p2) != 2:
        raise ValueError("Invalid vector points")
    if n_phi is not None:
        if not isinstance(n_phi, int) or n_phi < 3:
            raise ValueError("Invalid uniform angular partition")
    else:
        if not isinstance(phi, list) or len(phi) < 2:
            raise ValueError("Invalid custom angular partition")
        for i in range(1, len(phi)):
            if phi[i] <= phi[i - 1]:
                raise ValueError("Phi values should be increasing")
            if (phi[i] - phi[i - 1]) >= 180:
                raise ValueError("Phi step should be less then 180 degree")
        if phi[-1] - phi[0] > 360:
            raise ValueError("Phi values range should be <= 360 degree")
    # calculate phi's
    if n_phi is None:
        inp_phi = copy.deepcopy(phi)
    else:
        inp_phi = []
        for i in range(n_phi + 1):
            inp_phi.append(i * 360 / n_phi)

    c = com.grid3dcom.Revolve({"base": obj,
                               "p1": Point2(*p1),
                               "p2": Point2(*p2),
                               "phi": inp_phi,
                               "bt1": btype1,
                               "bt2": btype2,
                               "center_tri": not merge_central})
    try:
        flow.exec_command(c)
        return c._get_added_names()[2][0]
    except:
        raise ExecError("planar grid revolution")


def heal_grid(grid_id, simplify_boundary=30):
    """ Set of procedures for simplification of grid geometry

    Args:
       grid_id: identifier or (list of identifiers) of the grid

    Kwargs:
       simplify_boundary: angle (deg) in [0, 180].
          If this parameter is non-negative then edges which

          * are boundary,
          * belong to the same grid cell,
          * form an angle no more then ``simplify_boundary`` degree
          * elimination of intermediate vertices will not lead
            to cell degeneration

          will be merged. If ``simplify_boundary=0`` then only edges
          lying on the same line are considered, ``simplify_boundary=180``
          leads to merging of all doubled boundary edges,
          ``simplify_boundary=-1`` ignores this simplification option.

    Returns:
       None

    Raises:
       ValueError, hmscript.ExecError

    """
    if isinstance(grid_id, list):
        for g in grid_id:
            heal_grid(g, simplify_boundary)
        return
    sb = simplify_boundary
    if not isinstance(sb, (int, float, long)) or sb > 180 or sb < -1:
        raise ValueError("Invalid simplify_boundaries option: %s" % str(sb))

    c = com.gridcom.HealGrid({"name": grid_id, "simp_bnd": sb})
    try:
        flow.exec_command(c)
    except:
        raise ExecError('heal_grid')
