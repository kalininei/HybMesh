import ctypes as ct
from . import cport
import g2
from proc import (ccall, ccall_cb, list_to_c, concat, supplement,
                  CBoundaryNames)


def free_grid3(obj):
    ccall(cport.g3_free, obj)


def deepcopy(obj):
    ret = ct.c_void_p()
    ccall(cport.g3_deepcopy, obj, ct.byref(ret))
    return ret


def concatenate(objs):
    objs = list_to_c(objs, "void*")
    nobjs = ct.c_int(len(objs))
    ret = ct.c_void_p()
    ccall(cport.g3_concatenate, nobjs, objs, ct.byref(ret))
    return ret


def dims(obj):
    ret = (ct.c_int * 4)()
    ccall(cport.g3_dims, obj, ret)
    return list(ret)


def bnd_dims(obj):
    ret = (ct.c_int * 3)()
    ccall(cport.g3_bnd_dims, obj, ret)
    return list(ret)


def move(obj, dx, dy, dz):
    dx = (ct.c_double * 3)(dx, dy, dz)
    ccall(cport.g3_move, obj, dx)


def scale(obj, xpc, ypc, zpc, px, py, pz):
    p0 = (ct.c_double * 3)(px, py, pz)
    pc = (ct.c_double * 3)(xpc, ypc, zpc)
    ccall(cport.g3_scale, obj, pc, p0)


def raw_data(obj, what):
    ret = None
    d = dims(obj)
    if what == 'btypes':
        ret = (ct.c_int * d[2])()
        ccall(cport.g3_tab_btypes, obj, ret)
    else:
        raise Exception('unknown what: %s' % what)
    return ret


def point_by_index(obj, index):
    ret = (ct.c_double * 3)()
    ccall(cport.g3_point_at, obj, ct.c_int(index), ct.byref(ret))
    return list(ret)


def volume(obj):
    ret = ct.c_double()
    ccall(cport.g3_volume, obj, ct.byref(ret))
    return ret.value


def bnd_area(obj):
    ret = ct.c_double()
    ccall(cport.g3_bnd_area, obj, ret)
    return ret.value


def extract_surface(obj):
    ret = ct.c_void_p()
    ccall(cport.g3_extract_surface, obj, ct.byref(ret))
    return ret


def extrude(g2obj, zvals, bbot, btop, bside=None):
    d2 = g2.dims(g2obj)
    nzvals = ct.c_int(len(zvals))
    zvals = list_to_c(zvals, float)
    bbot = list_to_c(supplement(bbot, d2[2]), int)
    btop = list_to_c(supplement(btop, d2[2]), int)
    bside = ct.c_int(bside if bside is not None else -1)
    ret = ct.c_void_p()
    ccall(cport.g3_extrude, g2obj, nzvals, zvals, bbot, btop, bside,
          ct.byref(ret))
    return ret


def revolve(g2obj, vec, phivals, center_tri, bt1, bt2):
    vec = list_to_c(concat(vec), float)
    nphivals = ct.c_int(len(phivals))
    phivals = list_to_c(phivals, float)
    center_tri = ct.c_int(center_tri)
    bt1 = ct.c_int(bt1)
    bt2 = ct.c_int(bt2)
    ret = ct.c_void_p()
    ccall(cport.g3_revolve, g2obj, vec, nphivals, phivals,
          center_tri, bt1, bt2, ct.byref(ret))
    return ret


def tetrahedral_fill(sobjs, constrs, pts, pt_sizes, cb):
    nsobjs = ct.c_int(len(sobjs))
    sobjs = list_to_c(sobjs, 'void*')
    nconstrs = ct.c_int(len(constrs))
    constrs = list_to_c(constrs, 'void*')
    npts = ct.c_int(len(pts))
    pts = list_to_c(concat(pts), float)
    pt_sizes = list_to_c(pt_sizes, float)
    ret = ct.c_void_p()
    ccall_cb(cport.g3_tetrahedral_fill, cb,
             nsobjs, sobjs, nconstrs, constrs, npts, pts, pt_sizes,
             ct.byref(ret))
    return ret


def merge(obj1, obj2, cb=None):
    ret = ct.c_void_p()
    ccall_cb(cport.g3_merge, cb, obj1, obj2, ct.byref(ret))
    return ret


def to_msh(obj, fname, btypes, per_data, cb=None):
    """ per_data : [bnd_per-0, bnd_shadow-0,
                    pnt_per-0 as [x, y, z], pnt_shadow-0, ...]
    """
    fname = fname.encode('utf-8')
    btypes = CBoundaryNames(btypes)
    n_per_data = ct.c_int(len(per_data) / 4)
    it = iter(per_data)
    tmp = []
    while 1:
        try:
            tmp.append(next(it))
            tmp.append(next(it))
            tmp.extend(next(it))
            tmp.extend(next(it))
        except StopIteration:
            break
    per_data = list_to_c(tmp, float)
    ccall_cb(cport.g3_to_msh, cb, obj, fname, btypes, n_per_data, per_data)


def to_tecplot(obj, fname, btypes, cb=None):
    fname = fname.encode('utf-8')
    btypes = CBoundaryNames(btypes)
    ccall_cb(cport.g3_to_tecplot, cb, obj, fname, btypes)


def to_vtk(obj, fname, cb):
    fname = fname.encode('utf-8')
    ccall_cb(cport.g3_to_vtk, cb, obj, fname)


def to_gmsh(obj, fname, btypes, cb=None):
    fname = fname.encode('utf-8')
    btypes = CBoundaryNames(btypes)
    ccall_cb(cport.g3_to_gmsh, cb, obj, fname, btypes)


def surface_to_vtk(obj, fname, cb):
    fname = fname.encode('utf-8')
    ccall_cb(cport.g3_surface_to_vtk, cb, obj, fname)


def to_hm(doc, node, obj, name, fmt, afields, cb):
    name = name.encode('utf-8')
    naf = ct.c_int(len(afields))
    af = (ct.c_char_p * len(afields))()
    for i in range(len(afields)):
        af[i] = afields[i].encode('utf-8')
    ccall_cb(cport.g3_to_hm, cb, doc, node, obj, name, fmt, naf, af)
