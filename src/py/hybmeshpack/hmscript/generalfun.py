"general geometry and flow status functions"
from hybmeshpack.hmscript import flow
from hybmeshpack.hmscript import ExecError
from hybmeshpack import progdata
from hybmeshpack import com


def copy_geom(objs):
    """ Creates deep copies of geometry objects

    :param objs: identifier or list of identifiers of objects to copy

    :returns: list of identifiers of copied objects if objs is alist,
              single identifier of copied object otherwise
    """
    if isinstance(objs, list):
        return map(copy_geom, objs)

    try:
        c = com.objcom.CopyGeom({"names": [objs]})
        flow.exec_command(c)
        ret = c.added_contours2() + c.added_grids2() + c.added_surfaces3() +\
            c.added_grids3()
        return ret[0] if len(ret) > 0 else None
    except:
        raise ExecError("copy_geom")


def move_geom(objs, dx, dy, dz=0.):
    """ Moves a list of objects

    :param objs: identifier or list of identifiers of moving objects

    :param float dx:

    :param float dy:

    :param float dz: shifts in x, y and z direction. Z moves take place only
       for 3d objects

    :returns: None
    """
    ob = objs if isinstance(objs, list) else [objs]
    try:
        c = com.objcom.MoveGeom({"names": ob, "dx": dx, "dy": dy, "dz": dz})
        flow.exec_command(c)
    except:
        raise ExecError("move_geom")


def scale_geom(objs, xpc=100., ypc=100., zpc=100., refp=[0.0, 0.0, 0.0]):
    """ Scales objects

    :param objs: identifier or list of identifiers of scaling objects

    :param float xpc:

    :param float ypc:

    :param float zpc: percentages of scaling in x, y and z directions

    :param list-of-float refp: reference point as ``[x, y, z]`` which stays
       fixed after transformation

    :returns: None

    if **objs** contains only 2d objects **zpc** is ignored and could be
    ommitted, **refp** could by given as ``[x, y]``.
    """
    # for backward compatibility check if zpc was omitted
    if isinstance(zpc, list):
        refp = zpc
        zpc = 100.
    if len(refp) == 2:
        refp = [refp[0], refp[1], 0.]

    try:
        ob = objs if isinstance(objs, list) else [objs]
        c = com.objcom.ScaleGeom({"names": ob, "xpc": xpc, "ypc": ypc,
                                  "zpc": zpc, "p0": refp})
        flow.exec_command(c)
    except:
        raise ExecError("scale_geom")


def rotate_geom(objs, angle, pc=[0.0, 0.0]):
    """ Rotates group of 2d objects

    :param objs: identifier or list of identifiers of rotating 2d objects

    :param float angle: degree of rotation. Positive angle corresponds to
       counterclockwise rotation

    :param list-of-float pc: center of rotation

    :returns: None
    """
    try:
        ob = objs if isinstance(objs, list) else [objs]
        c = com.objcom.RotateGeom({"names": ob, "angle": angle,
                                   "p0": pc})
        flow.exec_command(c)
    except:
        raise ExecError("rotate_geom")


def reflect_geom(objs, pnt1, pnt2):
    """ Makes a reflection of 2d geometry objects over a  given line

    :param objs: identifier or list of identifiers of 2d objects to reflect

    :param list-of-float pnt1:

    :param list-of-float pnt2: points in [x, y] format which define a line
        to reflect over

    :returns: None
    """
    try:
        ob = objs if isinstance(objs, list) else [objs]
        c = com.objcom.ReflectGeom({"names": ob, "p1": pnt1, "p2": pnt2})
        flow.exec_command(c)
    except:
        raise ExecError("reflect_geom")


def remove_geom(objs):
    """ Completely removes object or list of objects

    :param objs: identifier or list of identifiers of removing objects

    :returns: None
    """
    try:
        ob = objs if isinstance(objs, list) else [objs]
        c = com.objcom.RemoveGeom({"names": ob})
        flow.exec_command(c)
    except:
        raise ExecError("remove_geom")


def remove_all():
    """ Completely removes all geometry objects and boundary types

    :returns: None
    """
    try:
        c = com.objcom.RemoveAll({})
        flow.exec_command(c)
    except:
        raise ExecError("remove_all")


def remove_all_but(objs):
    """ Removes all geometry objects except for listed ones

    :param objs: identifier or list of identifiers of objects
       which should NOT be removed

    :returns: None
    """
    if not isinstance(objs, list):
        objs = [objs]

    all_obj = flow.receiver.get_names()
    all_obj = [x for x in all_obj if x not in objs]
    remove_geom(all_obj)


def check_compatibility(vers, policy=1):
    """ Checks version compatibility. Notifies if current version
    of hymbesh is not fully compatible with input version.

    :param str vers: version in ``"0.1.2"`` format

    :param int policy: if versions are incompatible then:

       * 0 - do nothing (return False)
       * 1 - report warning to cout
       * 2 - raise Exception

    :returns: False if versions are incompatible, True otherwise.

    """
    versc = progdata.HybMeshVersion.current()
    versi = progdata.HybMeshVersion(vers)
    # last checked for 0.2.1 version
    ret = True
    if (versi > versc):
        ret = False
    if (versi < progdata.HybMeshVersion("0.2.1")):
        ret = False
    if not ret:
        message = "Current version of hybmesh %s is not " \
            "fully compatible with version %s" % (versc, vers)
        if policy == 1:
            print "WARNING:", message
        elif policy == 2:
            raise Exception(message)
    return ret


def registered_contours():
    """ Returns list of all contour identifiers
    """
    return flow.receiver.get_contour2_names()


def registered_grids():
    """ Returns list of all 2d grid identifiers
    """
    return flow.receiver.get_grid2_names()


def registered_grids3d():
    """ Returns list of all 3d grid identifiers
    """
    return flow.get_receiver().get_grid3_names()


def registered_surfaces():
    """ Returns list of all surface identifiers
    """
    return flow.receiver.get_surface3_names()


def registered_btypes():
    """ Returns list of all boundary types as list of (index, name) tuples
    """
    ret = []
    for k, v in flow.receiver.get_zone_types().iteritems():
        ret.append((k, v))
    return ret


def add_boundary_type(index, name="boundary1"):
    """ Register boundary type name.

    :param int index: index of boundary (>0)

    :param str name: user defined name of the boundary

    :returns: integer boundary identifier

    If boundary with ``index`` already exists it will be overwritten.
    Name of the boundary should be unique, if name already exists it will
    be changed automatically.

    """
    c = com.contcom.EditBoundaryType({"index": index, "name": name})
    flow.exec_command(c)
    return index
