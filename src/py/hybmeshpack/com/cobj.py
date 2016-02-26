' procedures for communication with c libraries'
import ctypes as ct
import os
from hybmeshpack.gdata.grid2 import Grid2
from hybmeshpack.gdata.contour2 import Contour2
from hybmeshpack.basic.geom import Point2
from hybmeshpack import progdata


def cport_lib():
    libfn = progdata.get_lib_fn('cport')
    return ct.cdll.LoadLibrary(os.path.abspath(libfn))


def grid_to_c(g):
    " return c pointer to a python grid"
    lib_fa = cport_lib()
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
    return lib_fa.grid_construct(ct.c_int(npt), ct.c_int(ncls), c_pnt, c_cls)


def free_c_grid(cgrid):
    lib_fa = cport_lib()
    lib_fa.grid_free(cgrid)


def cont_to_c(c):
    """ return c pointer to a contours2.Contour2 object
        works for contours of crossgrid library.
        Will be removed
    """
    lib_fa = cport_lib()
    #fill points
    c_pnt = (ct.c_double * (2 * c.n_points()))()
    for i, p in enumerate(c.points):
        c_pnt[2 * i] = p.x
        c_pnt[2 * i + 1] = p.y
    #fill edges
    edg = c.edges_points()
    c_edg = (ct.c_int * (2 * c.n_edges()))()
    for i, ce in enumerate(edg):
        c_edg[2 * i] = ce[0]
        c_edg[2 * i + 1] = ce[1]

    return lib_fa.contour_construct(
        ct.c_int(c.n_points()),
        ct.c_int(c.n_edges()), c_pnt, c_edg)


def grid_from_c(c_gr):
    "builds a grid from a c object"
    lib_fa = cport_lib()
    #c types
    ct_pint = ct.POINTER(ct.c_int)
    ct_ppint = ct.POINTER(ct.POINTER(ct.c_int))
    ct_pd = ct.POINTER(ct.c_double)
    ct_ppd = ct.POINTER(ct.POINTER(ct.c_double))
    #c data allocation
    npt, ned, ncl = ct.c_int(), ct.c_int(), ct.c_int()
    pt, ed, cdims, ced = ct_pd(), ct_pint(), ct_pint(), ct_pint()
    #call c function
    lib_fa.grid_get_edges_info.argtypes = [
        ct.c_void_p, ct_pint, ct_pint,
        ct_pint, ct_ppd, ct_ppint, ct_ppint, ct_ppint]
    lib_fa.grid_get_edges_info(
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
    lib_fa.grid_free_edges_info.argtypes = [ct_ppd, ct_ppint, ct_ppint,
                                            ct_ppint]
    lib_fa.grid_free_edges_info(ct.byref(pt), ct.byref(ed), ct.byref(cdims),
                                ct.byref(ced))

    return ret


def cont2_to_c(cont):
    """ return c pointer to a contours2.AbstractContour2 objects
    """
    lib_c2 = cport_lib()
    # arguments
    c_npnt = ct.c_int(cont.n_points())
    c_pnts = (ct.c_double * (2 * c_npnt.value))()
    c_nedges = ct.c_int(cont.n_edges())
    for i, p in enumerate(cont.points):
        c_pnts[2 * i] = p.x
        c_pnts[2 * i + 1] = p.y

    c_edges = (ct.c_int * (2 * c_nedges.value))()
    for i, [i0, i1, b] in enumerate(cont.edges_points()):
        c_edges[2 * i] = i0
        c_edges[2 * i + 1] = i1

    return lib_c2.create_ecollection_container(c_npnt, c_pnts,
                                               c_nedges, c_edges)


def cont2_from_c(c_cont):
    """ returns contours2.Contour2 from c pointer to ecollection """
    if c_cont == 0:
        return None
    lib_fa = cport_lib()
    c_eds, c_pts = ct.POINTER(ct.c_int)(), ct.POINTER(ct.c_double)()
    npts, neds = ct.c_int(), ct.c_int()
    lib_fa.ecollection_edges_info(c_cont,
                                  ct.byref(npts), ct.byref(neds),
                                  ct.byref(c_pts), ct.byref(c_eds))
    # from c arrays to python lists
    pts = []
    edcon = []
    for i in range(0, 2 * npts.value, 2):
        pts.append(Point2(c_pts[i], c_pts[i + 1]))
    for i in range(0, 2 * neds.value, 2):
        edcon.append([c_eds[i], c_eds[i + 1]])
    lib_fa.free_ecollection_edges_info(c_pts, c_eds)
    ret = Contour2.create_from_point_set(pts, edcon)
    return ret


def free_cont2(cont):
    """ frees pointer to ecollection container"""
    lib_fa = cport_lib()
    lib_fa.free_ecollection_container(cont)


def get_skewness(grid, threshold):
    """ reports skewness of the grid ->
           {'max_skew': float, 'max_skew_cell': int,
            'bad_cells': [cell indicies: int],
            'bad_skew': [cell skew: float]}
        returns {} if error
    """
    cgrid = grid_to_c(grid)
    lib_fa = cport_lib()

    ms = ct.c_double()
    msc = ct.c_int()
    num_bs = ct.c_int()
    bc = (ct.c_int * grid.n_cells())()
    bs = (ct.c_double * grid.n_cells())()
    cret = lib_fa.report_skewness(
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


def clip_domain(cc1, cc2, op, simplify):
    """ cc1, cc2 - c pointers to ecollection containers
    op is {1: union, 2: difference, 3: intersection, 4: xor}
    simplify (bool)
    Returns pointer to resulting container or None
    """
    lib_fa = cport_lib()
    c_op = ct.c_int(op)
    c_simp = ct.c_int(1 if simplify else 0)
    res = lib_fa.domain_clip(cc1, cc2, c_op, c_simp)
    if (res == 0):
        return None
    else:
        return res
