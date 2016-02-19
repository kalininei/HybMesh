from HybMeshPyPack import com
from HybMeshPyPack.hmscript import flow
from HybMeshPyPack.basic.geom import Point2


def remove_geom(objs):
    """ Completely removes object

    Args:
       obj: identifier of object or list of identifiers

    """
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.RemoveGeom({"names": ob})
    flow.exec_command(c)


def remove_all():
    """ Completely removes all grids, contours and btypes
    """
    flow.to_zero_state()


def move_geom(objs, dx, dy):
    """ Moves list of objects

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

       xpc, ypc (float): percentages of scaling in x, y directions

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
