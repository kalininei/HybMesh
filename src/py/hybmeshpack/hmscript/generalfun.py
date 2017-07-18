"general geometry and flow status functions"
import copy
from hybmeshpack.hmscript import flow, hmscriptfun
from hybmeshpack import progdata
from hybmeshpack import com
from datachecks import (icheck, ListOr1, UListOr1, Point2D, ZType, ACont2D,
                        OneOf, Float, String, AObject, APoint, NoneOr, List,
                        UInt, Func, InvalidArgument, Dict)


@hmscriptfun
def copy_geom(objs):
    """ Creates deep copies of geometry objects

    :param objs: identifier or list of identifiers of objects to copy

    :returns: list of identifiers of copied objects in oder prescribed by
       input list
    """
    icheck(0, ListOr1(AObject()))
    if not isinstance(objs, list):
        objs = [objs]

    c = com.objcom.CopyGeom({"names": objs})
    flow.exec_command(c)
    return c.odered_output()


@hmscriptfun
def move_geom(objs, dx, dy, dz=0.):
    """ Moves a list of objects

    :param objs: identifier or list of identifiers of moving objects

    :param float dx:

    :param float dy:

    :param float dz: shifts in x, y and z direction. Z moves take place only
       for 3d objects

    :returns: None
    """
    icheck(0, UListOr1(AObject()))
    icheck(1, Float())
    icheck(2, Float())
    icheck(3, Float())

    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.MoveGeom({"names": ob, "dx": dx, "dy": dy, "dz": dz})
    flow.exec_command(c)


@hmscriptfun
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
    icheck(0, UListOr1(AObject()))
    icheck(1, Float(grthan=0.))
    icheck(2, Float(grthan=0.))
    icheck(3, Float(grthan=0.))
    icheck(4, APoint())

    if len(refp) == 2:
        refp = [refp[0], refp[1], 0.]

    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.ScaleGeom({"names": ob, "xpc": xpc, "ypc": ypc,
                              "zpc": zpc, "p0": refp})
    flow.exec_command(c)


@hmscriptfun
def rotate_geom(objs, angle, pc=[0.0, 0.0]):
    """ Rotates group of 2d objects

    :param objs: identifier or list of identifiers of rotating 2d objects

    :param float angle: degree of rotation. Positive angle corresponds to
       counterclockwise rotation

    :param list-of-float pc: center of rotation

    :returns: None
    """
    icheck(0, UListOr1(ACont2D()))
    icheck(1, Float())
    icheck(2, Point2D())
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.RotateGeom({"names": ob, "angle": angle,
                               "p0": pc})
    flow.exec_command(c)


@hmscriptfun
def reflect_geom(objs, pnt1, pnt2):
    """ Makes a reflection of 2d geometry objects over a  given line

    :param objs: identifier or list of identifiers of 2d objects to reflect

    :param list-of-float pnt1:

    :param list-of-float pnt2: points in [x, y] format which define a line
        to reflect over

    :returns: None
    """
    icheck(0, UListOr1(ACont2D()))
    icheck(1, Point2D())
    icheck(2, Point2D(noteq=[pnt1]))

    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.ReflectGeom({"names": ob, "p1": pnt1, "p2": pnt2})
    flow.exec_command(c)


@hmscriptfun
def remove_geom(objs):
    """ Completely removes object or list of objects

    :param objs: identifier or list of identifiers of removing objects

    :returns: None
    """
    icheck(0, UListOr1(AObject()))
    ob = objs if isinstance(objs, list) else [objs]
    c = com.objcom.RemoveGeom({"names": ob})
    flow.exec_command(c)


@hmscriptfun
def remove_all():
    """ Completely removes all geometry objects and boundary types

    :returns: None
    """
    c = com.objcom.RemoveAll({})
    flow.exec_command(c)


@hmscriptfun
def remove_all_but(objs):
    """ Removes all geometry objects except for listed ones

    :param objs: identifier or list of identifiers of objects
       which should NOT be removed

    :returns: None
    """
    icheck(0, UListOr1(AObject()))

    if not isinstance(objs, list):
        objs = [objs]

    all_obj = flow.receiver.get_names()
    all_obj = [x for x in all_obj if x not in objs]
    remove_geom(all_obj)


@hmscriptfun
def check_compatibility(vers, policy=1):
    """ Checks version compatibility. Notifies if current version
    of hymbesh is not fully compatible with input version.

    :param str vers: version in ``"0.1.2"`` format

    :param int policy: if versions are incompatible then:

       * 0 - do nothing (return False)
       * 1 - report warning to cout
       * 2 - raise ExecError

    :returns: False if versions are incompatible, True otherwise.

    """
    icheck(0, String())
    icheck(1, OneOf(1, 2, 3))

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


@hmscriptfun
def registered_contours():
    """ Returns list of all contour identifiers
    """
    return flow.receiver.get_contour2_names()


@hmscriptfun
def registered_grids():
    """ Returns list of all 2d grid identifiers
    """
    return flow.receiver.get_grid2_names()


@hmscriptfun
def registered_grids3d():
    """ Returns list of all 3d grid identifiers
    """
    return flow.receiver.get_grid3_names()


@hmscriptfun
def registered_surfaces():
    """ Returns list of all surface identifiers
    """
    return flow.receiver.get_surface3_names()


@hmscriptfun
def registered_btypes():
    """ Returns list of all boundary types as list of (index, name) tuples
    """
    ret = []
    for k, v in flow.receiver.get_zone_types().iteritems():
        ret.append((k, v))
    return ret


@hmscriptfun
def add_boundary_type(index, name="boundary1"):
    """ Register boundary type name.

    :param int index: index of boundary (>0)

    :param str name: user defined name of the boundary

    :returns: integer boundary identifier

    If boundary with ``index`` already exists it will be overwritten.
    Name of the boundary should be unique, if name already exists it will
    be changed automatically.

    """
    icheck(0, ZType())
    icheck(1, String())

    c = com.contcom.EditBoundaryType({"index": index, "name": name})
    flow.exec_command(c)
    return index


def _setbt_b2grid(g, kv):
    iconts = g.raw_data('bnd')
    for k in kv.keys():
        for i, v in enumerate(kv[k]):
            kv[k][i] = iconts[v]


def _setbt_args_g2g3btps(g, btps, kv):
    _setbt_args_c2s3btps(g.contour(), btps, kv)
    _setbt_b2grid(g, kv)


def _setbt_args_g2bfun(g, bfun, kv):
    _setbt_args_c2bfun(g.contour(), bfun, kv)
    _setbt_b2grid(g, kv)


def _setbt_args_g3bfun(g, bfun, kv):
    _setbt_args_s3bfun(g.surface(), bfun, kv)
    _setbt_b2grid(g, kv)


def _setbt_args_c2s3btps(g, btps, kv):
    bold = g.raw_data('bt')
    for i, b1, b2 in zip(range(len(bold)), bold, btps):
        if b1 == b2:
            continue
        if b2 not in kv:
            kv[b2] = []
        kv[b2].append(i)


def _setbt_args_s3bfun(g, bfun, kv):
    bold = g.raw_data('bt')
    cnt = g.raw_data('face_center')
    it = iter(cnt)
    for i, b1, x, y, z in zip(range(len(bold)), bold, it, it, it):
        b2 = bfun(x, y, z, b1)
        if b2 is None or b2 == b1:
            continue
        if b2 not in kv:
            kv[b2] = []
        kv[b2].append(i)


def _setbt_args_c2bfun(g, bfun, kv):
    bold = g.raw_data('bt')
    vert = g.raw_data('vert')
    ev = g.raw_data('edge_vert')
    for i, b1 in enumerate(bold):
        e0, e1 = ev[2*i], ev[2*i+1]
        b2 = bfun(vert[2*e0], vert[2*e0+1],
                  vert[2*e1], vert[2*e1+1], b1)
        if b2 is None or b2 == b1:
            continue
        if b2 not in kv:
            kv[b2] = []
        kv[b2].append(i)


@hmscriptfun
def set_boundary_type(obj, btps=None, bfun=None, bdict=None):
    """ Mark geometrical object with boundary types.

    :param obj: geometric object identifier

    :param btps: single identifier for the whole object or list
       of identifiers for each boundary segment.

    :param  bfun: function which returns boundary type taking segment
       coordinates and old boundary type as arguments.

    :param bdict: {btype: [list-of-segment indicies]} dictionary
       which maps boundary type with object segments indicies

    Only one of **btps**, **bfun**, **bdict** arguments should be defined.

    **bfun** signature is:

       * ``(x0, y0, x1, y1, bt) -> btype`` for 2D objects, where
         *x0, y0, x1, y1* are edge end point coordinates,
         bt - old boundary type
       * ``(xc, yc, zc, bt) -> btype`` for 3D objects, where
         *xc, yc, zc* - approximate face center coordinates,
         bt - old boundary type

    If **obj** is a grid then only boundary segments will be passed
    to **bfun** function and **btps** list entries will
    be associated with boundary segments only.
    However **bdict** entries should contain global edge or face indicies.

    Example:

      .. literalinclude:: ../../testing/py/fromdoc/ex_setbtype.py
          :start-after: START OF EXAMPLE
          :end-before: END OF EXAMPLE

    """
    icheck(0, AObject())
    icheck(1, NoneOr(ListOr1(ZType())))
    icheck(3, NoneOr(Dict(ZType(), List(UInt()))))
    t = flow.receiver.whatis(obj)
    if t in ['g2', 'c2']:
        icheck(2, NoneOr(Func(narg=5)))
    else:
        icheck(2, NoneOr(Func(narg=4)))
    if [btps, bfun, bdict].count(None) != 2:
        raise InvalidArgument(
            "One of [btps, bfun, bdict] should be not None")

    args = {'name': obj, 'whole': None, 'btypes': {}}

    if isinstance(btps, int):
        args['whole'] = btps
    elif bdict is not None:
        args['btypes'] = bdict
    else:
        g = flow.receiver.get_object(obj)
        if t == 'g2':
            if btps is not None:
                _setbt_args_g2g3btps(g, btps, args['btypes'])
            if bfun is not None:
                _setbt_args_g2bfun(g, bfun, args['btypes'])
        if t == 'g3':
            if btps is not None:
                _setbt_args_g2g3btps(g, btps, args['btypes'])
            if bfun is not None:
                _setbt_args_g3bfun(g, bfun, args['btypes'])
        if t == 'c2':
            if btps is not None:
                _setbt_args_c2s3btps(g, btps, args['btypes'])
            if bfun is not None:
                _setbt_args_c2bfun(g, bfun, args['btypes'])
        if t == 's3':
            if btps is not None:
                _setbt_args_c2s3btps(g, btps, args['btypes'])
            if bfun is not None:
                _setbt_args_s3bfun(g, bfun, args['btypes'])
    c = com.objcom.SetBType(args)
    flow.exec_command(c)
