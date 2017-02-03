import ctypes as ct
from . import cport
from proc import (ccall, ccall_cb, list_to_c, free_cside_array, move_to_static,
                  CBoundaryNames, concat, supplement)


def dims(obj):
    ret = (ct.c_int * 3)()
    ccall(cport.g2_dims, obj, ret)
    return list(ret)


def area(obj):
    ret = ct.c_double(0)
    ccall(cport.g2_area, obj, ct.byref(ret))
    return ret.value


def bnd_dims(obj):
    ret = (ct.c_int * 2)()
    ccall(cport.g2_bnd_dims, obj, ret)
    return list(ret)


def bnd_length(obj):
    ret = ct.c_double()
    ccall(cport.g2_bnd_length, obj, ct.byref(ret))
    return ret.value


def skewness(obj, threshold):
    """ reports skewness of the grid ->
           {'max_skew': float, 'max_skew_cell': int,
            'bad_cells': [cell indicies: int],
            'bad_skew': [cell skew: float]}
        returns {} if error
    """
    maxskew = ct.c_double()
    maxskewindex = ct.c_int()
    badnum = ct.c_int()
    badindex = ct.POINTER(ct.c_int)()
    badvals = ct.POINTER(ct.c_double)()

    ccall(cport.g2_skewness, obj, threshold,
          ct.byref(maxskew), ct.byref(maxskewindex),
          ct.byref(badnum), ct.byref(badindex), ct.byref(badvals))

    ret = {}
    ret['max_skew'] = maxskew.value
    ret['max_skew_cell'] = maxskewindex.value
    ret['bad_cells'] = badindex[:badnum.value]
    ret['bad_skew'] = badvals[:badnum.value]

    free_cside_array(badindex, int)
    free_cside_array(badvals, float)
    return ret


def deepcopy(obj):
    ret = ct.c_void_p()
    ccall(cport.g2_deepcopy, obj, ct.byref(ret))
    return ret


def free_grid2(obj):
    ccall(cport.g2_free, obj)


def concatenate(objs):
    objs = list_to_c(objs, "void*")
    nobjs = ct.c_int(len(objs))
    ret = ct.c_void_p()
    ccall(cport.g2_concatenate, nobjs, objs, ct.byref(ret))
    return ret


def move(obj, dx, dy):
    dx = (ct.c_double * 2)(dx, dy)
    ccall(cport.g2_move, obj, dx)


def scale(obj, xpc, ypc, px, py):
    p0 = (ct.c_double * 2)(px, py)
    pc = (ct.c_double * 2)(xpc, ypc)
    ccall(cport.g2_scale, obj, pc, p0)


def reflect(obj, x0, y0, x1, y1):
    v0 = (ct.c_double * 2)(x0, y0)
    v1 = (ct.c_double * 2)(x1, y1)
    ccall(cport.g2_reflect, obj, v0, v1)


def rotate(obj, x0, y0, angle):
    p0 = (ct.c_double * 2)(x0, y0)
    a = ct.c_double(angle)
    ccall(cport.g2_rotate, obj, p0, a)


def set_bnd(obj, bndlist, for_all_edges):
    nb = dims(obj)[1] if for_all_edges else bnd_dims(obj)[1]
    bndlist = list_to_c(supplement(bndlist, nb), int)
    for_all_edges = ct.c_int(for_all_edges)
    ccall(cport.g2_set_btypes, obj, bndlist, for_all_edges)


def extract_contour(obj):
    ret = ct.c_void_p()
    ccall(cport.g2_extract_contour, obj, ct.byref(ret))
    return ret


def raw_data(obj, what):
    ret = None
    d = dims(obj)
    if what == 'btypes':
        ret = (ct.c_int * d[1])()
        ccall(cport.g2_tab_btypes, obj, ret)
    elif what == 'vertices':
        ret = ((ct.c_double * 2) * d[0])()
        ccall(cport.g2_tab_vertices, obj, ret)
    elif what == 'edge-vert':
        ret = ((ct.c_int * 2) * d[1])()
        ccall(cport.g2_tab_edgevert, obj, ret)
    elif what == 'cellsizes':
        ret = (ct.c_int * d[2])()
        ccall(cport.g2_tab_cellsizes, obj, ret)
    elif what == 'cell-vert':
        nret, ret2 = ct.c_int(), ct.POINTER(ct.c_int)()
        ccall(cport.g2_tab_cellvert, obj, ct.byref(nret), ct.byref(ret2))
        ret = move_to_static(nret, ret2, int)
    elif what == 'cell-edge':
        nret, ret2 = ct.c_int(), ct.POINTER(ct.c_int)()
        ccall(cport.g2_tab_celledge, obj, ct.byref(nret), ct.byref(ret2))
        ret = move_to_static(nret, ret2, int)
    elif what == 'centers':
        ret = ((ct.c_double * d[2]) * 2)()
        ccall(cport.g2_tab_centers, obj, ret)
    elif what == 'bedges':
        nret, ret2 = ct.c_int(), ct.POINTER(ct.c_int)()
        ccall(cport.g2_tab_bedges, obj, ct.byref(nret), ct.byref(ret2))
        ret = move_to_static(nret, ret2, int)
    else:
        raise ValueError('unknown what: %s' % what)
    return ret


def point_by_index(obj, index):
    ret = ct.c_int()
    ccall(cport.g2_point_at, obj, ct.c_int(index), ct.byref(ret))
    return ret.value


def closest_points(obj, pts, proj):
    npts = ct.c_int(len(pts))
    pts = list_to_c(concat(pts), 'float')
    ret = (ct.c_double * len(pts))()
    if proj == "vertex":
        proj = ct.c_int(0)
    elif proj == "edge":
        proj = ct.c_int(1)
    else:
        raise ValueError
    ccall(cport.g2_closest_points, obj, npts, pts, proj, ret)
    it = iter(ret)
    return [[a, b] for a, b in zip(it, it)]


def from_points_edges(points, eds):
    """ eds = [[p0, p1, cleft, cright, btype], ....] """
    npoints = ct.c_int(len(points))
    points = list_to_c(concat(points), float)
    neds = ct.c_int(len(eds))
    eds = list_to_c(concat(eds), int)
    ret = ct.c_void_p()

    ccall(cport.g2_from_points_edges, npoints, points, neds, eds,
          ct.byref(ret))

    return ret


def from_points_cells(points, cls, bedges):
    """
    cls = [[p1, p2, p3, ...], [p1, p2, ....]]
    bedges: [[p1, p2, btype], ...]
    """
    npoints = ct.c_int(len(points))
    points = list_to_c(concat(points), float)
    ncells = ct.c_int(len(cls))
    cellsizes = list_to_c(map(len, cls), int)
    cellvert = list_to_c(concat(cls), int)
    nbedges = ct.c_int(len(bedges))
    bedges = list_to_c(concat(bedges), int)
    ret = ct.c_void_p()

    ccall(cport.g2_from_points_cells, npoints, points,
          ncells, cellsizes, cellvert,
          nbedges, bedges, ct.byref(ret))

    return ret


def build_rect_grid(xdata, ydata, bnds):
    nx = ct.c_int(len(xdata))
    xdata = list_to_c(xdata, float)
    ny = ct.c_int(len(ydata))
    ydata = list_to_c(ydata, float)
    bnds = list_to_c(supplement(bnds, 4), 'int')
    ret = ct.c_void_p()
    ccall(cport.g2_rect_grid, nx, xdata, ny, ydata, bnds, ct.byref(ret))
    return ret


def build_circ_grid(p0, rdata, adata, istrian, bnd):
    p0 = list_to_c(p0, float)
    nr = ct.c_int(len(rdata))
    rdata = list_to_c(rdata, float)
    na = ct.c_int(len(adata))
    adata = list_to_c(adata, float)
    istrian = ct.c_int(istrian)
    bnd = ct.c_int(bnd)
    ret = ct.c_void_p()
    ccall(cport.g2_circ_grid, p0, nr, rdata, na, adata, istrian, bnd,
          ct.byref(ret))
    return ret


def build_ring_grid(p0, rdata, adata, bnds):
    p0 = list_to_c(p0, float)
    nr = ct.c_int(len(rdata))
    rdata = list_to_c(rdata, float)
    na = ct.c_int(len(adata))
    adata = list_to_c(adata, float)
    bnds = list_to_c(supplement(bnds, 2), int)
    ret = ct.c_void_p()
    ccall(cport.g2_ring_grid, p0, nr, rdata, na, adata, bnds, ct.byref(ret))
    return ret


def build_tri_grid(verts, nedge, bnds):
    verts = list_to_c(concat(verts), float)
    nedge = ct.c_int(nedge)
    bnds = list_to_c(supplement(bnds, 3), int)
    ret = ct.c_void_p()
    ccall(cport.g2_tri_grid, verts, nedge, bnds, ct.byref(ret))
    return ret


def regular_hex_grid(area, crad, strict):
    area = list_to_c(area, float)
    crad = ct.c_double(crad)
    areatype = "hex" if len(area) == 3 else "rect"
    strict = ct.c_int(strict)
    ret = ct.c_void_p()
    ccall(cport.g2_hex_grid, areatype, area, crad, strict, ct.byref(ret))
    return ret


def unstructed_fill(domain, constraint, embpts, filler):
    nembpts = ct.c_int(len(embpts))
    embpts = list_to_c(concat(embpts), float)
    ret = ct.c_void_p()
    ccall(cport.g2_unstructured_fill, domain, constraint,
          nembpts, embpts, filler, ct.byref(ret))
    return ret


def custom_rectangular_grid(algo, left, bottom, right, top,
                            herw, rinvalid, cb):
    herw = list_to_c(supplement(herw, 4), float)
    rinvalid = ct.c_int(rinvalid)
    ret = ct.c_void_p()
    ccall_cb(cport.g2_custom_rect_grid, cb,
             algo, left, bottom, right, top,
             herw, rinvalid, ct.byref(ret))
    return ret


def circ4grid(algo, p0, rad, step, sqrside, rcoef):
    p0 = list_to_c(p0, float)
    rad = ct.c_double(rad)
    step = ct.c_double(step)
    sqrside = ct.c_double(sqrside)
    rcoef = ct.c_double(rcoef)
    ret = ct.c_void_p()
    ccall(cport.g2_circ4grid, algo, p0, rad, step, sqrside, rcoef,
          ct.byref(ret))
    return ret


def stripe_grid(obj, partition, tipalgo, bnd, cb):
    npartition = ct.c_int(len(partition))
    partition = list_to_c(partition, float)
    bnd = list_to_c(supplement(bnd, 4), int)
    ret = ct.c_void_p()
    ccall_cb(cport.g2_stripe_grid, cb, obj,
             npartition, partition, tipalgo, bnd, ct.byref(ret))
    return ret


def map_grid(base_obj, target_obj, base_points, target_points,
             snap, bt_from_contour, algo, is_reversed, rinvalid, cb):
    npoints = min(len(base_points), len(target_points))
    base_points = base_points[:npoints]
    target_points = target_points[:npoints]
    npoints = ct.c_int(npoints)
    base_points = list_to_c(concat(base_points), float)
    target_points = list_to_c(concat(target_points), float)
    bt_from_contour = ct.c_int(bt_from_contour)
    is_reversed = ct.c_int(is_reversed)
    rinvalid = ct.c_int(rinvalid)
    ret = ct.c_void_p()

    ccall_cb(cport.g2_map_grid, cb, base_obj, target_obj,
             npoints, base_points, target_points, snap,
             bt_from_contour, algo, is_reversed, rinvalid, ct.byref(ret))
    return ret


def simplify_grid_boundary(obj, angle):
    angle = ct.c_double(angle)
    ret = ct.c_void_p()
    ccall(cport.g2_simplify_bnd, obj, angle, ct.byref(ret))
    return ret


def convex_cells(obj, angle):
    angle = ct.c_double(angle)
    ret = ct.c_void_p()
    ccall(cport.g2_convex_cells, obj, angle, ct.byref(ret))
    return ret


def grid_excl_cont(obj, cont, isinner, cb):
    isinner = ct.c_int(isinner)
    ret = ct.c_void_p()
    ccall_cb(cport.g2_exclude_cont, cb, obj, cont, isinner, ct.byref(ret))
    return ret


def unite_grids(obj1, obj2, buf, fixbnd, emptyholes, angle0, filler, cb):
    buf = ct.c_double(buf)
    fixbnd = ct.c_int(fixbnd)
    emptyholes = ct.c_int(emptyholes)
    angle0 = ct.c_double(angle0)
    ret = ct.c_void_p()
    ccall_cb(cport.g2_unite_grids, cb,
             obj1, obj2, buf, fixbnd, emptyholes, angle0, filler,
             ct.byref(ret))
    return ret


def to_msh(obj, fname, btypes, per_data, cb=None):
    fname = fname.encode('utf-8')
    btypes = CBoundaryNames(btypes)
    n_per_data = ct.c_int(len(per_data) / 3)
    per_data = list_to_c(per_data, int)
    ccall(cport.g2_to_msh, obj, fname, btypes, n_per_data, per_data)


def to_tecplot(obj, fname, btypes, cb=None):
    fname = fname.encode('utf-8')
    btypes = CBoundaryNames(btypes)
    ccall(cport.g2_to_tecplot, obj, fname, btypes)


def to_hm(doc, node, obj, name, fmt, afields, cb=None):
    name = name.encode('utf-8')
    naf = ct.c_int(len(afields))
    af = (ct.c_char_p * len(afields))()
    for i in range(len(afields)):
        af[i] = afields[i].encode('utf-8')
    ccall(cport.g2_to_hm, doc, node, obj, name, fmt, naf, af)


def boundary_layer_grid(opt, cb=None):
    """ ->grid2.
        opt - options dictionary object.
              Same as gridcom.BuildBoundaryGrid options['opt'] dictionary
              except for 'source' field is a proper Contour2.cdata object.
              All fields must be filled
        cb - callback of CB_CANCEL2 type
    """
    # prepare c input data
    class COptStruct(ct.Structure):
        def __init__(self, opt_entry):
            self.cont = opt_entry['source']
            n = len(opt_entry['partition'])
            self.Npartition = ct.c_int(n)
            self.partition = (ct.c_double * n)(*opt_entry['partition'])
            self.direction = ct.c_char_p(
                {1: 'LEFT', -1: 'RIGHT'}[opt_entry['direction']])
            self.mesh_cont = ct.c_char_p({
                0: 'NO', 1: 'KEEP_ORIGIN', 2: 'KEEP_SHAPE', 3: 'IGNORE_ALL',
                4: "INCREMENTAL",
            }[opt_entry['mesh_cont']])
            self.mesh_cont_step = ct.c_double(opt_entry['mesh_cont_step'])
            self.start = (ct.c_double * 2)(opt_entry['start'].x,
                                           opt_entry['start'].y)
            self.end = (ct.c_double * 2)(opt_entry['end'].x,
                                         opt_entry['end'].y)
            self.force_conformal = ct.c_int(
                1 if opt_entry['force_conf'] else 0)
            self.angle_range = (ct.c_double * 4)(
                opt_entry['algo_acute'],
                opt_entry['algo_right'],
                opt_entry['algo_straight'],
                opt_entry['algo_reentr'])
            self.step_start = ct.c_double(opt_entry['step_start'])
            self.step_end = ct.c_double(opt_entry['step_end'])

        _fields_ = [
            ('cont', ct.c_void_p),
            ('Npartition', ct.c_int),
            ('partition', ct.POINTER(ct.c_double)),
            ('direction', ct.c_char_p),
            ('mesh_cont', ct.c_char_p),
            ('mesh_cont_step', ct.c_double),
            ('start', ct.c_double * 2),
            ('end', ct.c_double * 2),
            ('force_conformal', ct.c_int),
            ('angle_range', ct.c_double * 4),
            ('step_start', ct.c_double),
            ('step_end', ct.c_double),
        ]

    # 2) build array of c structures
    nopt = ct.c_int(len(opt))
    c_opt_type = COptStruct * len(opt)
    c_opt = c_opt_type()
    for i, co in enumerate(opt):
        c_opt[i] = COptStruct(co)

    ret = ct.c_void_p()
    ccall_cb(cport.g2_boundary_layer, cb, nopt, c_opt, ct.byref(ret))
    return ret
