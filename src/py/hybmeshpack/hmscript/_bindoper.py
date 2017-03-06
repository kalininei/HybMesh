from hybmeshpack.hmscript import flow, hmscriptfun
from datachecks import (icheck, AObject, ZType, NoneOr, List, UInt, OneOf)
import c2oper
import o2info
import o3info
import generalfun
import ctypes as ct
import struct


# special functions for high-level interfaces
@hmscriptfun
def dims(gid):
    icheck(0, AObject())
    ob = flow.receiver.get_object(gid)
    return ob.dims()


@hmscriptfun
def partition_segment(start, end, hstart, hend, hinternal=[]):
    ret1 = c2oper.partition_segment(start, end, hstart, hend, hinternal)
    return (ct.c_double * len(ret1))(*ret1)


@hmscriptfun
def skewness(gid, threshold):
    ret1 = o2info.skewness(gid, threshold)
    total = len(ret1['bad_cells'])
    buf = ((ct.c_char * 12) * total)()
    for i, (c, s) in enumerate(zip(ret1['bad_cells'], ret1['bad_skew'])):
        struct.pack_into("=id", buf, 12*i, c, s)
    return buf


@hmscriptfun
def set_btypes(obj, bt, segm=None):
    icheck(0, AObject())
    icheck(1, ZType())
    icheck(2, NoneOr(List(UInt())))
    if segm is None:
        generalfun.set_boundary_type(obj, bt)
    else:
        generalfun.set_boundary_type(obj, bdict={bt: segm})


@hmscriptfun
def raw_data(obj, what):
    icheck(0, AObject())
    t = flow.receiver.whatis(obj)
    if t == 'c2':
        return o2info.tab_cont2(obj, what)
    elif t == 'g2':
        return o2info.tab_grid2(obj, what)
    elif t == 's3':
        return o3info.tab_surf3(obj, what)
    elif t == 'g3':
        return o3info.tab_grid3(obj, what)

@hmscriptfun
def stdout_verbosity(verb):
    icheck(0, OneOf(0, 1, 2, 3))
    flow.set_interface(flow.interface.reset(verb))
