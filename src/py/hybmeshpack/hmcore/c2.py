import ctypes as ct
from . import cport
from proc import (ccall, list_to_c, free_cside_array,
                  supplement, concat)


def deepcopy(obj):
    ret = ct.c_void_p()
    ccall(cport.c2_deepcopy, obj, ct.byref(ret))
    return ret


def free_cont2(obj):
    ccall(cport.c2_free, obj)


def area(obj):
    ret = ct.c_double()
    ccall(cport.c2_area, obj, ct.byref(ret))
    return ret.value


def dims(obj):
    ret = (ct.c_int * 2)()
    ccall(cport.c2_dims, obj, ret)
    return list(ret)


def length(obj):
    ret = ct.c_double()
    ccall(cport.c2_length, obj, ct.byref(ret))
    return ret.value


def move(obj, dx, dy):
    dx = (ct.c_double * 2)(dx, dy)
    ccall(cport.c2_move, obj, dx)


def scale(obj, xpc, ypc, px, py):
    p0 = (ct.c_double * 2)(px, py)
    pc = (ct.c_double * 2)(xpc, ypc)
    ccall(cport.c2_scale, obj, pc, p0)


def reflect(obj, x0, y0, x1, y1):
    v0 = (ct.c_double * 2)(x0, y0)
    v1 = (ct.c_double * 2)(x1, y1)
    ccall(cport.c2_reflect, obj, v0, v1)


def rotate(obj, x0, y0, angle):
    p0 = (ct.c_double * 2)(x0, y0)
    a = ct.c_double(angle)
    ccall(cport.c2_rotate, obj, p0, a)


def build_from_points(pts, force_closed, bnds):
    bnds = list_to_c(supplement(bnds, len(pts)), int)
    force_closed = ct.c_int(force_closed)
    npts = ct.c_int(len(pts))
    pts = list_to_c(concat(pts), float)
    ret = ct.c_void_p()
    ccall(cport.c2_frompoints, npts, pts, bnds, force_closed, ct.byref(ret))
    return ret


def simplify(obj, angle, keep_btypes):
    angle = ct.c_double(angle)
    keep_btypes = ct.c_int(keep_btypes)
    ccall(cport.c2_simplify_self, obj, angle, keep_btypes)


def quick_separate(obj):
    ret_nconts = ct.c_int(0)
    ret_conts = ct.POINTER(ct.c_void_p)()

    ccall(cport.c2_quick_separate,
          obj, ct.byref(ret_nconts), ct.byref(ret_conts))
    ret = []
    for i in range(ret_nconts.value):
        ret.append(ret_conts[i])
    free_cside_array(ret_conts, 'void*')
    return ret


def unite_contours(objs):
    nconts = ct.c_int(len(objs))
    objs = list_to_c(objs, 'void*')
    ret_cont = ct.c_void_p(0)
    ccall(cport.c2_unite, nconts, objs, ct.byref(ret_cont))
    return ret_cont


def decompose_contour(obj):
    ret = ct.POINTER(ct.c_void_p)()
    nret = ct.c_void_p()
    ccall(cport.c2_decompose, obj, ct.byref(nret), ct.byref(ret))
    a = [ret[i] for i in range(nret.value)]
    free_cside_array(ret, 'void*')
    return a


def spline(pts, bnds, nedges):
    bnds = list_to_c(supplement(bnds, len(pts)), int)
    npts = ct.c_int(len(pts))
    pts = list_to_c(concat(pts), float)
    nedges = ct.c_int(nedges)
    ret = ct.c_void_p()
    ccall(cport.c2_spline, npts, pts, nedges, bnds, ct.byref(ret))
    return ret


def clip_domain(obj1, obj2, op, simp):
    "op is [union, difference, intersection, xor]"
    simp = ct.c_int(simp)
    ret = ct.c_void_p()
    ccall(cport.c2_clip_domain, obj1, obj2, op, simp, ct.byref(ret))
    return ret


def contour_partition(obj, step, algo, a0, keepbnd,
                      nedges, crosses, keeppts, start, end):
    """ step: [double]*n, where n dependes on algo:
            algo = "const" => n = 1,
            algo = "ref_points" => n = 3*k
            algo = "ref_weights/lengths" => n = 2*k
        nedges: None or forced number of resulting edges
    """
    # options
    nstep = ct.c_int(len(step))
    step = list_to_c(step, "float")
    a0 = ct.c_double(a0)
    keepbnd = ct.c_int(keepbnd)
    nedges = ct.c_int(nedges) if nedges is not None else ct.c_int(-1)
    ncrosses = ct.c_int(len(crosses))
    crosses = list_to_c(crosses, "void*")
    nkeeppts = ct.c_int(len(keeppts))
    keeppts = list_to_c(concat(keeppts), float)
    start = list_to_c(start, float) if start is not None else None
    end = list_to_c(end, float) if end is not None else None

    # call
    ret = ct.c_void_p()
    ccall(cport.c2_partition, obj, algo,
          nstep, step, a0, keepbnd, nedges,
          ncrosses, crosses,
          nkeeppts, keeppts,
          start, end,
          ct.byref(ret))
    return ret


def matched_partition(obj, conts, pts, step, infdist, power, a0):
    nconts = ct.c_int(len(conts))
    conts = list_to_c(conts, "void*")
    npts = ct.c_int(len(pts) / 3)
    pts = list_to_c(pts, float)
    step = ct.c_double(step)
    infdist = ct.c_double(infdist)
    power = ct.c_double(power)
    a0 = ct.c_double(a0)

    ret = ct.c_void_p()
    ccall(cport.c2_matched_partition, obj, nconts, conts,
          npts, pts, step, infdist, power, a0,
          ct.byref(ret))
    return ret


def segment_partition(start, end, hstart, hend, hinternal):
    start = ct.c_double(start)
    end = ct.c_double(end)
    hstart = ct.c_double(hstart)
    hend = ct.c_double(hend)
    ninternal = ct.c_int(len(hinternal) / 2)
    hinternal = list_to_c(hinternal, float)
    nret = ct.c_int()
    ret = ct.POINTER(ct.c_double)()

    ccall(cport.c2_segment_partition, start, end, hstart, hend,
          ninternal, hinternal, ct.byref(nret), ct.byref(ret))
    r = ret[:nret.value]
    free_cside_array(ret, float)
    return r


def extract_subcontours(obj, plist):
    nplist = ct.c_int(len(plist))
    plist = list_to_c(concat(plist), float)
    ret = (ct.c_void_p * (nplist.value - 1))()
    ccall(cport.c2_extract_subcontours, obj, nplist, plist, ret)
    return list(ret)


def connect_subcontours(objs, fx, close, shift):
    nobjs = ct.c_int(len(objs))
    objs = list_to_c(objs, "void*")
    nfx = ct.c_int(len(fx))
    fx = list_to_c(fx, int)
    shift = ct.c_int(shift)
    ret = ct.c_void_p()
    ccall(cport.c2_connect_subcontours, nobjs, objs, nfx, fx, shift, close,
          ct.byref(ret))
    return ret


def set_bnd(obj, bnd):
    bnd = list_to_c(bnd, "int")
    ccall(cport.c2_set_btypes, obj, bnd)


def contour_type(obj):
    ret = ct.c_int()
    ccall(cport.c2_contour_type, obj, ct.byref(ret))
    ret = ret.value
    if ret == 0:
        return 'open'
    elif ret == 1:
        return 'closed'
    elif ret == 2:
        return 'mdomain'
    else:
        return 'compound'


def concatenate(objs):
    objs = list_to_c(objs, "void*")
    nobjs = ct.c_int(len(objs))
    ret = ct.c_void_p()
    ccall(cport.c2_concatenate, nobjs, objs, ct.byref(ret))
    return ret


def to_hm(doc, node, obj, name, fmt, cb):
    ccall(cport.c2_to_hm, doc, node, obj, name, fmt)


def closest_points(obj, pts, proj):
    npts = ct.c_int(len(pts))
    pts = list_to_c(concat(pts), 'float')
    ret = (ct.c_double * len(pts))()
    ccall(cport.c2_closest_points, obj, npts, pts, proj, ret)
    it = iter(ret)
    return [[a, b] for a, b in zip(it, it)]


def raw_data(obj, what):
    ret = None
    d = dims(obj)
    if what == 'btypes':
        ret = (ct.c_int * d[1])()
        ccall(cport.c2_tab_btypes, obj, ret)
    elif what == 'vertices':
        ret = ((ct.c_double * 2) * d[0])()
        ccall(cport.c2_tab_vertices, obj, ret)
    elif what == 'edge-vert':
        ret = ((ct.c_int * 2) * d[1])()
        ccall(cport.c2_tab_edgevert, obj, ret)
    else:
        raise Exception('unknown what: %s' % what)
    return ret


def point_by_index(obj, index):
    ret = (ct.c_double * 2)()
    ccall(cport.c2_point_at, obj, ct.c_int(index), ct.byref(ret))
    return list(ret)


def end_points(obj):
    start = (ct.c_double * 2)()
    end = (ct.c_double * 2)()
    ccall(cport.c2_end_points, obj, start, end)
    return list(start), list(end)
