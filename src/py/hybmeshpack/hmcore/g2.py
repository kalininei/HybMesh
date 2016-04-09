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
