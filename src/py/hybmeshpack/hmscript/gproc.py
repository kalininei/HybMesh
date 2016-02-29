from hybmeshpack import com
from hybmeshpack.hmscript import flow
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


def move_geom(objs, dx, dy):
    """ Moves a list of objects

    Args:
       objs: identifier or list of identifiers of moving objects

       dx, dy (float): shifts in x and y direction

    """
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.MoveGeom({"names": ob, "dy": dy, "dx": dx})
    flow.exec_command(c)


def rotate_geom(objs, angle, pc=[0.0, 0.0]):
    """ Rotates group of objects

    Args:
       objs: identifier or list of identifiers of rotating objects

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
       objs: identifier or list of identifiers of scaling objects

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
    elif isinstance(objs, str):
        c = com.objcom.CopyGeom({"names": [objs], "newnames": [objs]})
        flow.exec_command(c)
        if len(c._get_added_names()[0]) > 0:
            return c._get_added_names()[0]
        else:
            return c._get_added_names()[1]


def heal_grid(grid_id, simplify_boundary=30):
    """ Set of procedures for simplification of grid geometry

    Args:
       grid_id: identifier of the grid

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
    sb = simplify_boundary
    if not isinstance(sb, (int, float, long)) or sb > 180 or sb < -1:
        raise ValueError("Invalid simplify_boundaries option: %s" % str(sb))

    c = com.gridcom.HealGrid({"name": grid_id, "simp_bnd": sb})
    try:
        flow.exec_command(c)
    except:
        raise ExecError('heal_grid')
