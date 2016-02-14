from HybMeshPyPack import com, basic
import HybMeshPyPack.com.objcom
from HybMeshPyPack.hmscript import flow
from HybMeshPyPack.basic.geom import Point2


def RemoveGeom(obj):
    """ delete grid or contour:
        obj - string identifier of object
    """
    c = com.objcom.RemoveGeom({"names": [obj]})
    flow.exec_command(c)


def MoveGeom(objs, dx, dy):
    """ moves list of grids or contours
        objs - list of string identifiers of moving objects
        dx, dy - shifts in x and y direction
    """
    c = com.objcom.MoveGeom({"names": objs, "dy": dy, "dx": dx})
    flow.exec_command(c)


def RotateGeom(objs, angle, pc=[0.0, 0.0]):
    """ rotates group of grids and contours
       objs - list of string identifiers of moving objects
       angle - degree of rotation. Positive angle corresponds to
               counterclockwise rotation
       pc - center of rotation as [x, y]
    """
    c = com.objcom.RotateGeom({"names": objs, "angle": angle,
       "p0": Point2(*pc)})
    flow.exec_command(c)


def ScaleGeom(objs, xpc, ypc, refp=[0.0, 0.0]):
    """ scales grids and contours
       objs - list of string identifiers of scaled objects
       xpc, ypc - percentages of scalint in x, y directions
       pc - refference point as [x, y]
    """
    c = com.objcom.ScaleGeom({"names": objs, "xpc": xpc, "ypc": ypc,
        "p0": Point2(*refp)})
    flow.exec_command(c)


def CopyGeom(objs):
    """ creates deep copies of contours and grids
        objs - list of string identifiers of copied objects
        Returns list of identifiers of copied objects
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
