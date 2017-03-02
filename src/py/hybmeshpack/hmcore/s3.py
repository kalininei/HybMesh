import ctypes as ct
from . import cport
from proc import (ccall, ccall_cb, free_cside_array, list_to_c,
                  move_to_static, BndTypesDifference)


def deepcopy(obj):
    ret = ct.c_void_p()
    ccall(cport.s3_deepcopy, obj, ct.byref(ret))
    return ret


def free_surf3(obj):
    ccall(cport.s3_free, obj)


def area(obj):
    ret = ct.c_double()
    ccall(cport.s3_area, obj, ct.byref(ret))
    return ret.value


def volume(obj):
    ret = ct.c_double()
    ccall(cport.s3_volume, obj, ct.byref(ret))
    return ret.value


def dims(obj):
    ret = (ct.c_int * 3)()
    ccall(cport.s3_dims, obj, ret)
    return list(ret)


def move(obj, dx, dy, dz):
    dx = (ct.c_double * 3)(dx, dy, dz)
    ccall(cport.s3_move, obj, dx)


def scale(obj, xpc, ypc, zpc, px, py, pz):
    p0 = (ct.c_double * 3)(px, py, pz)
    pc = (ct.c_double * 3)(xpc, ypc, zpc)
    ccall(cport.s3_scale, obj, pc, p0)


def quick_separate(obj):
    nret = ct.c_int(0)
    ret = ct.POINTER(ct.c_void_p)()

    ccall(cport.s3_quick_separate, obj, ct.byref(nret), ct.byref(ret))
    ret2 = []
    for i in range(nret.value):
        ret2.append(ret[i])
    free_cside_array(ret, 'void*')
    return ret2


def concatenate(objs):
    objs = list_to_c(objs, "void*")
    nobjs = ct.c_int(len(objs))
    ret = ct.c_void_p()
    ccall(cport.s3_concatenate, nobjs, objs, ct.byref(ret))
    return ret


def raw_data(obj, what):
    ret = None
    d = dims(obj)
    if what == 'vert':
        ret = (ct.c_double * (3*d[0]))()
        ccall(cport.s3_tab_vertices, obj, ret)
    elif what == 'edge_vert':
        ret = (ct.c_int * (2*d[1]))()
        ccall(cport.s3_tab_edgevert, obj, ret)
    elif what == 'face_dim':
        ret = (ct.c_int * d[2])()
        ccall(cport.s3_tab_facedim, obj, ret)
    elif what == 'face_edge':
        nret, cret = ct.c_int(), ct.POINTER(ct.c_int)()
        ccall(cport.s3_tab_faceedge, obj, ct.byref(nret), ct.byref(cret))
        ret = move_to_static(nret.value, cret, int)
    elif what == 'face_vert':
        nret, cret = ct.c_int(), ct.POINTER(ct.c_int)()
        ccall(cport.s3_tab_facevert, obj, ct.byref(nret), ct.byref(cret))
        ret = move_to_static(nret.value, cret, int)
    elif what == 'bt':
        ret = (ct.c_int*d[2])()
        ccall(cport.s3_tab_btypes, obj, ret)
    elif what == 'center':
        ret = (ct.c_double*(3*d[2]))()
        ccall(cport.s3_tab_centers, obj, ret)
    else:
        raise Exception('unknown what: %s' % repr(what))
    return ret


def point_by_index(obj, index):
    ret = (ct.c_double * 3)()
    ccall(cport.s3_point_at, obj, ct.c_int(index), ct.byref(ret))
    return list(ret)


def assign_boundary_types(obj, bt):
    dataout = ct.POINTER(ct.c_int)()
    ccall(cport.s3_assign_boundary_types, obj, bt.data, ct.byref(dataout))
    return BndTypesDifference.from_cdata(dataout)


def to_hm(doc, node, obj, name, fmt, cb):
    name = name.encode('utf-8')
    ccall_cb(cport.s3_to_hm, cb, doc, node, obj, name, fmt)
