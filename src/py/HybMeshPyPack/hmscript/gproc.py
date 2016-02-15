from HybMeshPyPack import com, basic
import HybMeshPyPack.com.objcom
from HybMeshPyPack.hmscript import flow
from HybMeshPyPack.basic.geom import Point2


def RemoveGeom(obj):
    """ Completely removes object

    Args:
       obj: identifier of object

    """
    c = com.objcom.RemoveGeom({"names": [obj]})
    flow.exec_command(c)


def MoveGeom(objs, dx, dy):
    """ Moves list of objects

    Args:
       objs: list of identifiers of moving objects

       dx, dy (float): shifts in x and y direction

    """
    c = com.objcom.MoveGeom({"names": objs, "dy": dy, "dx": dx})
    flow.exec_command(c)


def RotateGeom(objs, angle, pc=[0.0, 0.0]):
    """ Rotates group of objects
    
    Args:
      objs: list of string identifiers of moving objects
      
      angle (float): degree of rotation. Positive angle corresponds to
      counterclockwise rotation

      pc (list-of-float): center of rotation
    """
    c = com.objcom.RotateGeom({"names": objs, "angle": angle,
       "p0": Point2(*pc)})
    flow.exec_command(c)


def ScaleGeom(objs, xpc, ypc, refp=[0.0, 0.0]):
    """ Scales objects

    Args:
       objs: list of string identifiers of scaled objects

       xpc, ypc (float): percentages of scaling in x, y directions

       refp (list-of-float): reference point which stays
       fixed after transformation
    """
    c = com.objcom.ScaleGeom({"names": objs, "xpc": xpc, "ypc": ypc,
        "p0": Point2(*refp)})
    flow.exec_command(c)


def CopyGeom(objs):
    """ Creates deep copies of objects

    Args:
       objs: list of identifiers of objects to copy

    Returns:
       list of identifiers of copied objects
    """
    if isinstance(objs, list):
        ret = []
        for s in objs:
            ret.append(CopyGeom(s))
        return ret
    elif isinstance(objs, str):
        c = com.objcom.CopyGeom({"names": [objs], "newnames": [objs]})
        flow.exec_command(c)
        if len(c._get_added_names()[0]) > 0:
            return c._get_added_names()[0][0]
        else:
            return c._get_added_names()[1][0]
