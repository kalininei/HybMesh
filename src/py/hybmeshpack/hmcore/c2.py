import ctypes as ct
from . import libhmcport, list_to_c
from hybmeshpack.basic.geom import Point2
from hybmeshpack.gdata.contour2 import Contour2


def cont2_to_c(cont):
    """ return c pointer to a contours2.AbstractContour2 objects
    """
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

    return libhmcport.create_ecollection_container(c_npnt, c_pnts,
                                                   c_nedges, c_edges)


def cont2_from_c(c_cont, c_bnd=None):
    """ returns contours2.Contour2 from c pointer to ecollection """
    if c_cont == 0:
        return None
    c_eds, c_pts = ct.POINTER(ct.c_int)(), ct.POINTER(ct.c_double)()
    npts, neds = ct.c_int(), ct.c_int()
    libhmcport.ecollection_edges_info(c_cont,
                                      ct.byref(npts), ct.byref(neds),
                                      ct.byref(c_pts), ct.byref(c_eds))
    # from c arrays to python lists
    pts = []
    edcon = []
    for i in range(0, 2 * npts.value, 2):
        pts.append(Point2(c_pts[i], c_pts[i + 1]))
    for i in range(0, 2 * neds.value, 2):
        edcon.append([c_eds[i], c_eds[i + 1]])
    libhmcport.free_ecollection_edges_info(c_pts, c_eds)
    ret = Contour2.create_from_point_set(pts, edcon)
    if c_bnd is not None:
        bmap = {}
        for i, b in enumerate(c_bnd):
            if b != 0:
                bmap[i] = b
        ret.set_edge_bnd(bmap)
    return ret


def free_cont2(cont):
    """ frees pointer to ecollection container"""
    libhmcport.free_ecollection_container(cont)


def clip_domain(cc1, cc2, op, simplify):
    """ cc1, cc2 - c pointers to ecollection containers
    op is {1: union, 2: difference, 3: intersection, 4: xor}
    simplify (bool)
    Returns pointer to resulting container or None
    """
    c_op = ct.c_int(op)
    c_simp = ct.c_int(1 if simplify else 0)
    res = libhmcport.domain_clip(cc1, cc2, c_op, c_simp)
    if (res == 0):
        return None
    else:
        return res


def contour_partition(c_cont, c_bt, c_step, algo, a0, keepbnd, nedges):
    """ c_step: ct.c_double() * n, where n dependes on algo:
    algo = "const" => n = 1,
    algo = "ref_points" => n = 3*k
    """
    # algo treatment
    if algo == "const":
        c_algo = ct.c_int(0)
        c_n = ct.c_int(0)
    elif algo == "ref_points":
        c_algo = ct.c_int(1)
        c_n = ct.c_int(len(c_step) / 3)
    else:
        raise ValueError("invalid partition algo")

    # prepare arrays for boundary output
    c_bnd = ct.POINTER(ct.c_int)()
    n_bnd = ct.c_int(0)
    c_ne = ct.c_int(-1 if nedges is None else nedges)

    ret = libhmcport.contour_partition(c_cont, c_bt, c_algo,
                                       c_n, c_step,
                                       ct.c_double(a0),
                                       ct.c_int(1 if keepbnd else 0), c_ne,
                                       ct.byref(n_bnd), ct.byref(c_bnd)
                                       )
    # copy c-side allocated bnd array to python side allocated
    c_bndret = (ct.c_int * (n_bnd.value))()
    for i in range(len(c_bndret)):
        c_bndret[i] = c_bnd[i]
    # clear c-side allocated array
    libhmcport.free_int_array(c_bnd)

    if ret == 0:
        raise Exception("contour partition failed")
    return ret, c_bndret


def spline(c_pnts, c_bt, nedges):
    # algo treatment
    c_bnd = ct.POINTER(ct.c_int)()
    n_bnd = ct.c_int(0)
    c_ne = ct.c_int(nedges)
    npnt = len(c_pnts) / 2
    nbt = len(c_bt)

    ret = libhmcport.spline(npnt, c_pnts, nbt, c_bt, c_ne,
                            ct.byref(n_bnd), ct.byref(c_bnd))
    # copy c-side allocated bnd array to python side allocated
    c_bndret = (ct.c_int * (n_bnd.value))()
    for i in range(len(c_bndret)):
        c_bndret[i] = c_bnd[i]
    # clear c-side allocated array
    libhmcport.free_int_array(c_bnd)

    if ret == 0 or ret is None:
        raise Exception("spline builder failed")
    return ret, c_bndret


def matched_partition(c_cont, c_src, c_pts, step, dist, pw, a0):
    c_step = ct.c_double(step)
    c_a0 = ct.c_double(a0)
    c_dist = ct.c_double(dist)
    c_nc = ct.c_int(len(c_src))
    c_src2 = list_to_c(c_src, 'void*')
    c_pw = ct.c_double(pw)
    c_np = ct.c_int(len(c_pts) / 3) if c_pts is not None else ct.c_int(0)

    ret = libhmcport.matched_partition(c_cont, c_nc, c_src2, c_np, c_pts,
                                       c_step, c_dist, c_pw, c_a0)
    if ret == 0:
        raise Exception("contour partition failed")
    return ret
