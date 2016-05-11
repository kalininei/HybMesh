import ctypes as ct
from . import libhmcport, list_to_c
from hybmeshpack.basic.geom import Point2 as Point2
from hybmeshpack.gdata.grid2 import Grid2 as Grid2


def grid_to_c(g):
    " return c pointer to a python grid"
    npt = g.n_points()
    ncls = g.n_cells()
    #fill points
    c_pnt = (ct.c_double * (2 * npt))()
    for i, p in enumerate(g.points):
        c_pnt[2 * i] = p.x
        c_pnt[2 * i + 1] = p.y
    #fill cells
    cls_a = g.cells_nodes_connect()
    cls_a_num = len(cls_a)
    for ce in cls_a:
        cls_a_num += len(ce)
    c_cls = (ct.c_int * cls_a_num)()
    ind = 0
    for ce in cls_a:
        c_cls[ind] = len(ce)
        for i in range(1, len(ce) + 1):
            c_cls[ind + i] = ce[i - 1]
        ind += (len(ce) + 1)
    return libhmcport.grid_construct(
        ct.c_int(npt), ct.c_int(ncls), c_pnt, c_cls)


def free_c_grid(cgrid):
    libhmcport.grid_free(cgrid)


def grid_from_c(c_gr):
    "builds a grid from a c object"
    #c types
    ct_pint = ct.POINTER(ct.c_int)
    ct_ppint = ct.POINTER(ct.POINTER(ct.c_int))
    ct_pd = ct.POINTER(ct.c_double)
    ct_ppd = ct.POINTER(ct.POINTER(ct.c_double))
    #c data allocation
    npt, ned, ncl = ct.c_int(), ct.c_int(), ct.c_int()
    pt, ed, cdims, ced = ct_pd(), ct_pint(), ct_pint(), ct_pint()
    #call c function
    libhmcport.grid_get_edges_info.argtypes = [
        ct.c_void_p, ct_pint, ct_pint,
        ct_pint, ct_ppd, ct_ppint, ct_ppint, ct_ppint]
    libhmcport.grid_get_edges_info(
        c_gr, ct.byref(npt), ct.byref(ned),
        ct.byref(ncl), ct.byref(pt), ct.byref(ed),
        ct.byref(cdims), ct.byref(ced))
    # ---- construct grid
    ret = Grid2()
    #points
    it = iter(pt)
    for i in range(npt.value):
        ret.points.append(Point2(next(it), next(it)))
    #edges
    it = iter(ed)
    for i in range(ned.value):
        ret.edges.append([next(it), next(it)])
    #cells
    it1, it2 = iter(cdims), iter(ced)
    for i in range(ncl.value):
        ied = [next(it2) for j in range(next(it1))]
        ret.cells.append(ied)
    #free c memory
    libhmcport.grid_free_edges_info.argtypes = [ct_ppd, ct_ppint, ct_ppint,
                                                ct_ppint]
    libhmcport.grid_free_edges_info(ct.byref(pt), ct.byref(ed),
                                    ct.byref(cdims), ct.byref(ced))
    return ret


def get_skewness(grid, threshold):
    """ reports skewness of the grid ->
           {'max_skew': float, 'max_skew_cell': int,
            'bad_cells': [cell indicies: int],
            'bad_skew': [cell skew: float]}
        returns {} if error
    """
    cgrid = grid_to_c(grid)

    ms = ct.c_double()
    msc = ct.c_int()
    num_bs = ct.c_int()
    bc = (ct.c_int * grid.n_cells())()
    bs = (ct.c_double * grid.n_cells())()
    cret = libhmcport.report_skewness(
        cgrid, ct.c_double(threshold),
        ct.byref(ms), ct.byref(msc), ct.byref(num_bs),
        bc, bs)
    free_c_grid(cgrid)
    ret = {}
    if (cret != 0):
        return ret

    ret['max_skew'] = ms.value
    ret['max_skew_cell'] = msc.value
    ret['bad_cells'] = []
    ret['bad_skew'] = []

    for i in range(num_bs.value):
        ret['bad_cells'].append(bc[i])
        ret['bad_skew'].append(bs[i])

    return ret


def boundary_types_to_c(grid):
    """ returns pointer to c-structure which holdes boundary types
    for input grid. This structure is only actual for current grid
    state and should be freed after usage
    """
    nd1, nd2, bt = [], [], []
    for cont in grid.boundary_contours():
        for e in cont:
            ebnd = grid.get_edge_bnd(e)
            nd1.append(grid.edges[e][0])
            nd2.append(grid.edges[e][1])
            bt.append(ebnd)
    c1 = list_to_c(nd1, int)
    c2 = list_to_c(nd2, int)
    c3 = list_to_c(bt, int)
    n = ct.c_int(len(nd1))
    return libhmcport.set_grid2_boundary_types(n, c1, c2, c3)


def free_boundary_types(bt):
    """ frees memory allocated by boundary_types_to_c procedure """
    libhmcport.free_grid2_boundary_types(bt)


def custom_rectangular_grid(algo, c_left, c_bot, c_right, c_top, cb):
    """ algo: one of ['linear', 'inverse_laplace',
                      'direct_laplace', 'orthogonal']
        c_*: c allocated contours.
        returns c allocated 2d grid or raises
        After procedure points coordinates of c_* data may be changed.
    """
    if algo == 'linear':
        c_algo = ct.c_int(0)
    elif algo == 'inverse_laplace':
        c_algo = ct.c_int(1)
    elif algo == 'direct_laplace':
        c_algo = ct.c_int(2)
    elif algo == 'orthogonal':
        c_algo = ct.c_int(3)
    else:
        raise ValueError("Invalid custom rectangular grid algo")

    args = (c_algo, c_left, c_bot, c_right, c_top)
    cb.initialize(libhmcport.custom_rectangular_grid, args)
    cb.execute_command()
    res = cb.get_result()
    if res == 0:
        raise Exception('Custom rectangular grid builder failed')
    return res


def to_msh(c_g, fname, c_btypes, c_bnames, c_periodic):
    c_fname = fname.encode('utf-8')
    if c_periodic is not None:
        n_periodic = ct.c_int(len(c_periodic) / 3)
    else:
        n_periodic = ct.c_int(0)

    args = (c_g, c_fname, c_btypes, c_bnames, n_periodic, c_periodic)
    res = libhmcport.export_msh_grid(*args)

    if res != 0:
        raise Exception("msh 2d grid export failed")


def to_tecplot(c_g, fname, c_btypes, c_bnames):
    c_fname = fname.encode('utf-8')
    args = (c_g, c_fname, c_btypes, c_bnames)
    res = libhmcport.export_tecplot_grid(*args)
    if res != 0:
        raise Exception("tecplot 2d grid export failed")


def unite_grids(c_g1, c_g2, buf, fix_bnd, empty_holes, an0, cb):
    """ adds g2 to g1. Returns new c-grid.
        cb -- Callback.CB_CANCEL2 callback object

        raises Exception if failed
    """
    c_buf = ct.c_double(buf)
    c_an0 = ct.c_double(an0)
    c_fix = ct.c_int(1) if fix_bnd else ct.c_int(0)
    c_eh = ct.c_int(1) if empty_holes else ct.c_int(0)
    args = (c_g1, c_g2, c_buf, c_fix, c_eh, c_an0)

    libhmcport.cross_grids_wcb.restype = ct.c_void_p
    cb.initialize(libhmcport.cross_grids_wcb, args)
    cb.execute_command()
    ret = cb.get_result()

    #if result was obtained (no errors, no cancel)
    if ret != 0:
        return ret
    else:
        raise Exception("unite_grids failed")


def grid_excl_cont(c_grd, c_cnt, is_inner, cb):
    """ ->c_grid
        Returns a c-grid with excluded contour area (inner or outer)
        raises if fails
    """
    c_isinner = ct.c_int(1 if is_inner else 0)
    args = (c_grd, c_cnt, c_isinner)

    cb.initialize(libhmcport.grid_exclude_cont_wcb, args)
    cb.execute_command()
    res = cb.get_result()
    if res != 0:
        return res
    else:
        raise Exception("Grid exclusion failed")


def circ4grid(algo, c_p0, rad, step, sqrside, rcoef):
    """ -> c_grid or raise """
    if algo == "linear":
        c_algo = ct.c_int(0)
    elif algo == "laplace":
        c_algo = ct.c_int(1)
    elif algo == "orthogonal_circ":
        c_algo = ct.c_int(2)
    elif algo == "orthogonal_rect":
        c_algo = ct.c_int(3)
    else:
        raise ValueError("Unknown algorithm")

    res = libhmcport.circ4grid(
        c_algo, c_p0, ct.c_double(rad),
        ct.c_double(step), ct.c_double(sqrside),
        ct.c_double(rcoef))

    if res != 0:
        return res
    else:
        raise Exception("Quadrangular grid in circular area "
                        "grid builder failed")


def map_grid(c_grid, c_cont, c_gpoints, c_cpoints, snap, algo,
             return_invalid, cb):
    """maps grid on cont using gpoints, cpoints as basis.
       snap = "no", "add_vertices", "shift_vertices"
       algo = "inverse_laplace", "direct_laplace"
       bt = "from_grid", "from_contour"
    """
    n = ct.c_int(len(c_gpoints) / 2)
    # snap
    s = ct.c_int(0)
    if snap == "add_vertices":
        s = ct.c_int(2)
    elif snap == "shift_vertices":
        s = ct.c_int(3)
    elif snap == "no":
        s = ct.c_int(1)
    # algo
    a = ct.c_int(0)
    if algo == "direct_laplace":
        a = ct.c_int(1)
    elif algo == "inverse_laplace":
        a = ct.c_int(2)
    #invalid
    inv = ct.c_int(1 if return_invalid else 0)

    args = (c_grid, c_cont, n, c_gpoints, c_cpoints, s, a, inv)
    cb.initialize(libhmcport.build_grid_mapping, args)
    cb.execute_command()
    cret = cb.get_result()
    if cret == 0:
        raise Exception("Grid mapping failed")
    else:
        return cret
